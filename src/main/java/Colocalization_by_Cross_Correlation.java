


import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.*;

import net.imagej.axis.*;
import net.imagej.display.ImageDisplayService;
import net.imagej.ops.OpService;

import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.math.ImgMath;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.outofbounds.OutOfBoundsFactory;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;
import net.imglib2.util.Intervals;
import net.imglib2.util.Util;
import net.imglib2.view.Views;

import org.apache.commons.io.FileUtils;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;

import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.io.IOService;
import org.scijava.log.LogService;
import org.scijava.plot.*;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.table.Table;
import org.scijava.table.Tables;
import org.scijava.ui.DialogPrompt;
import org.scijava.ui.UIService;
import org.scijava.ui.swing.viewer.plot.jfreechart.XYPlotConverter;
import org.scijava.util.*;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.util.*;
import java.util.stream.Collectors;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, menuPath = "Analyze>Colocalization>Colocalization by Correlation")
public class Colocalization_by_Cross_Correlation implements Command{

    @Parameter
    private LogService logService;

    @Parameter
    private StatusService statusService;

    @Parameter
    private UIService uiService;

    @Parameter
    private DatasetService datasetService;

    @Parameter
    private ImageDisplayService imageDisplayService;

    @Parameter
    private PlotService plotService;

    @Parameter
    private IOService ioService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private OpService ops;

    @Parameter(label = "Image 1: ", description = "This is the image which will be randomized during Costes randomization", persist = false)
    private Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    private Dataset dataset2;

    @Parameter(label = "No mask (not recommended)?", description = "When checked, performs Costes randomization over the entire image, regardless of what image is selected below.", callback = "maskCallback")
    private boolean maskAbsent;

    @Parameter(label = "Mask: ", description = "The mask over which pixels of image 1 will be randomized. This is important, more details at: imagej.github.io/Colocalization_by_Cross_Correlation", required = false, persist = false)
    private Dataset maskDataset;

    @Parameter(label = "Cycle count: ", description = "The number of Costes randomization cycles to perform. Recommend at least 3, more for sparse signal.", min = "1")
    private long cycles;

    @Parameter(label = "Significant digits: ")
    private int significantDigits;

    @Parameter(label = "Show intermediate images? ", description = "Shows images of numerous steps throughout the algorithm. More details at: imagej.github.io/Colocalization_by_Cross_Correlation")
    private boolean showIntermediates;

    @Parameter(label = "Output directory (leave blank for none):", description = "The directory to automatically save all generated output, including the intermediate images if the \"Show Intermediates\" box is checked", required = false, style="directory")
    private File saveFolder;

    @Parameter(type = ItemIO.OUTPUT)
    private Dataset ContributionOf1, ContributionOf2, timeCorrelationHeatMap;

    @Parameter(type = ItemIO.OUTPUT, label = "Table of correlations")
    private Table correlationTable;

    @Parameter(type = ItemIO.OUTPUT, label = "CC Results")
    private Table results;

    @Parameter(type = ItemIO.OUTPUT)
    private XYPlot plot;

    @Parameter(type = ItemIO.OUTPUT, label = "Summary of Results")
    private String notes;

    private Dataset [] intermediates;

    private double sigDigits;

    private String [] intermediateNames = {"Original CC result", "Costes randomized image", "Subtracted CC result", "Gaussian-modified CC result"};

    private double [] scale;

    private String statusBase = "";

    public Colocalization_by_Cross_Correlation() {
    }

    @Override
    public void run(){

        //region Error checking
        if(dataset1.numDimensions() != dataset2.numDimensions() || dataset1.getHeight() != dataset2.getHeight() || dataset1.getWidth() != dataset2.getWidth() || dataset1.getDepth() != dataset2.getDepth() || dataset1.getFrames() != dataset2.getFrames()){
            logService.error("All image dimensions (XYZ, and time) must match");
            return;
        }
        if(!maskAbsent && (dataset1.numDimensions() != maskDataset.numDimensions() || dataset1.getHeight() != maskDataset.getHeight() || dataset1.getWidth() != maskDataset.getWidth() || dataset1.getDepth() != maskDataset.getDepth() || dataset1.getFrames() != maskDataset.getFrames())){
            logService.error("All mask dimensions (XYZ, and time) must match image dimensions");
            return;
        }
        if(dataset1.getChannels() > 1 || dataset2.getChannels() > 1 || (!maskAbsent && maskDataset.getChannels() > 1)){
            logService.error("Multi-channel images are not supported, requires separate channels");
            return;
        }
        //endregion

        statusService.showStatus("Initializing plugin data");

        //region Plugin initialization (mostly creating datasets)
        RadialProfiler radialProfile;
        sigDigits = Math.pow(10.0, significantDigits);

        if(maskAbsent){
            maskDataset = dataset1.duplicateBlank();
            LoopBuilder.setImages(maskDataset).multiThreaded().forEachPixel(SetOne::setOne);
        }

        SCIFIOConfig config = new SCIFIOConfig();
        config.writerSetFailIfOverwriting(false);

        // Cannot use duplicateBlank() for creating the upcoming images, as they need to be 32-bit Float images
        CalibratedAxis [] calibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
        AxisType [] axisTypes = new AxisType[dataset1.numDimensions()];
        for (int i = 0; i < dataset1.numDimensions(); ++i) {
            calibratedAxes[i] = dataset1.axis(i);
            axisTypes[i] = dataset1.axis(i).type();
        }

        ContributionOf1 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset1.getName(), axisTypes);
        ContributionOf1.setAxes(calibratedAxes);
        ContributionOf2 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset2.getName(), axisTypes);
        ContributionOf2.setAxes(calibratedAxes);


        if(showIntermediates) {
            intermediates = new Dataset[4];
            for (int i = 0; i < 4; i++) {
               intermediates[i] = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), intermediateNames[i] + " of " + dataset1.getName(), axisTypes);
                intermediates[i].setAxes(calibratedAxes);
            }
        }
        //endregion

        //region Single frame analysis
        if(dataset1.getFrames() == 1) {
            scale = new double[dataset1.numDimensions()];
            for (int i = 0; i < scale.length; i++) {
                scale[i] = dataset1.averageScale(i);
            }

            try {
                radialProfile = new RadialProfiler(dataset1, scale);
            } catch (Exception e) {
                e.printStackTrace();
                return;
            }
            try{colocalizationAnalysis(dataset1.duplicate(), dataset2.duplicate(), maskDataset, radialProfile, ContributionOf1, ContributionOf2, intermediates);}
            catch (Exception e){
                e.printStackTrace();
                throw e;
            }

            if(showIntermediates){
                for (int i = 0; i < 4; i++) {
                    uiService.show(intermediates[i]);
                }
            }


            SeriesStyle oCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("blue"), LineStyle.SOLID, MarkerStyle.NONE);
            SeriesStyle sCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("green"), LineStyle.SOLID, MarkerStyle.NONE);
            SeriesStyle gaussStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("red"), LineStyle.DASH, MarkerStyle.NONE);

            plot = plotService.newXYPlot();

            XYSeries gaussData = plot.addXYSeries();
            gaussData.setValues(radialProfile.gaussCurveMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).boxed().collect(Collectors.toList()), new ArrayList<>(radialProfile.gaussCurveMap.values()));
            gaussData.setStyle(gaussStyle);
            gaussData.setLabel("Gaussian Fit");

            XYSeries sCorrPlotData = plot.addXYSeries();
            sCorrPlotData.setValues(radialProfile.sCorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).boxed().collect(Collectors.toList()), new ArrayList<>(radialProfile.sCorrMap.values()));
            sCorrPlotData.setStyle(sCorrStyle);
            sCorrPlotData.setLabel("Subtracted CC");

            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(radialProfile.oCorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).boxed().collect(Collectors.toList()), new ArrayList<>(radialProfile.oCorrMap.values()));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");

            plot.xAxis().setLabel("Distance (" + getUnitType() + ")");

            plot.yAxis().setLabel("Cross correlation");

            plot.xAxis().setManualRange(0, Math.max(radialProfile.gaussFit[1] + (5 * radialProfile.gaussFit[2]), 0.01));

            plot.setTitle("Correlation of images");

            if(radialProfile.confidence < 15){
                notes = (notes == null ? "" : notes) + "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.\n\n\n";
            }

            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Gaussian Height");
            rowHeaders.add("Confidence");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(getSigDigits(radialProfile.gaussFit[1]));
            resultsList.add(getSigDigits(radialProfile.gaussFit[2]));
            resultsList.add(getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFit[1]))));
            resultsList.add(getSigDigits(radialProfile.confidence));

            results = Tables.wrap(resultsList, "", rowHeaders);

            List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            RadialProfiler finalRadialProfile = radialProfile;
            radialProfile.oCorrMap.keySet().stream().forEachOrdered((d) -> {
                LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
                row.put("Distance (" + getUnitType() +")", (getSigDigits(d.doubleValue())));
                row.put("Original CC", getSigDigits(finalRadialProfile.oCorrMap.get(d)));
                row.put("Subtracted CC", getSigDigits(finalRadialProfile.sCorrMap.get(d)));
                row.put("Gaussian fit", getSigDigits(finalRadialProfile.gaussCurveMap.get(d)));
                correlationTableList.add(row);
            });
            correlationTable = Tables.wrap(correlationTableList, null);

            if(saveFolder != null){
                if(!saveFolder.exists() || !saveFolder.canWrite()){
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);
                    datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + "\\" + ContributionOf1.getName() + ".tif", config);
                    datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + "\\" + ContributionOf2.getName() + ".tif", config);

                    File plotout = new File(saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".png");
                    XYPlotConverter converter = new XYPlotConverter();
                    ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth()*2, plot.getPreferredHeight()*2);

                    ioService.save(results,saveFolder.getAbsolutePath() + "\\" + "CC Results.csv" );
                    ioService.save(correlationTable, saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".csv");

                    String summary = (notes == null ? "" : notes) + "Fit a gaussian curve to the cross-correlation of: \n\""
                            + dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\"\n using the mask \n\"" +
                            (maskAbsent? "No mask selected" : maskDataset.getName()) +
                            "\":\n\nMean (" + getUnitType() +"): " + getSigDigits(radialProfile.gaussFit[1]) +
                            "\nStandard deviation: " + getSigDigits(radialProfile.gaussFit[2]) +
                            "\nGaussian height:" + getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFit[1]))) +
                            "\nConfidence: " + getSigDigits(radialProfile.confidence);

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Summary.txt"), summary, (Charset) null);

                    if(showIntermediates){
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + "\\" + intermediate.getName() + ".tif", config);
                        }
                    }

                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion

        //region Multi-frame analysis
        else{
            long[] min = new long[dataset1.numDimensions()];
            long[] max = dataset1.dimensionsAsLongArray();
            int timeAxis = dataset1.dimensionIndex(Axes.TIME);

            Optional<CalibratedAxis> calibratedTime = dataset1.axis(Axes.TIME);

            scale = new double[dataset1.numDimensions()-1];
            int j = 0;
            for (int i = 0; i < dataset1.numDimensions(); i++) {
                if(i != timeAxis) {
                    scale[j++] = dataset1.averageScale(i);
                }
            }

            for (int i = 0; i < max.length; i++) {
                max[i] = max[i]-1;
            }

            min[timeAxis] = 0;
            max[timeAxis] = 0;

            RandomAccess<RealType<?>> correlationAccessor = null;

            ArrayList<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            ArrayList<String> correlationTablesRowNames = new ArrayList<>();

            double highestConfidence = 0;
            double highestConMean = 0;
            double highestConSD = 0;
            long highestConFrame = 0;
            double highestCCvalue = 0;


            RandomAccessibleInterval [] intermediatesViewsPasser = null;
            if(showIntermediates){
                intermediatesViewsPasser = new RandomAccessibleInterval[4];
            }

            //Duplicate datasets to not modify originals when applying masks later
            Dataset dataset1copy = dataset1.duplicate();
            Dataset dataset2copy = dataset2.duplicate();

            Dataset tempHeatMap = null;

            for (long i = 0; i < dataset1.getFrames(); i++) {
                statusService.showProgress((int)i, (int)dataset1.getFrames());
                //StatusService message update is performed in colocalizationAnalysis method
                statusBase = "Frame " + (i+1) + " - ";
                min[timeAxis] = i;
                max[timeAxis] = i;

                RandomAccessibleInterval temp1 = Views.dropSingletonDimensions(Views.interval(dataset1copy, min, max));
                RandomAccessibleInterval temp2 = Views.dropSingletonDimensions(Views.interval(dataset2copy, min, max));
                RandomAccessibleInterval masktemp = Views.dropSingletonDimensions(Views.interval(maskDataset, min, max));

                if(showIntermediates){
                    for (int m = 0; m < 4; m++) {
                        intermediatesViewsPasser[m] = Views.dropSingletonDimensions(Views.interval(intermediates[m], min, max));
                    }
                }

                try {
                    radialProfile = new RadialProfiler(temp1, scale);
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }

                try{colocalizationAnalysis(datasetService.create(temp1), datasetService.create(temp2), datasetService.create(masktemp), radialProfile, Views.dropSingletonDimensions(Views.interval(ContributionOf1, min, max)), Views.dropSingletonDimensions(Views.interval(ContributionOf2, min, max)), intermediatesViewsPasser);}
                catch (Exception e){
                    e.printStackTrace();
                    throw e;
                }

                if(tempHeatMap == null) {
                    tempHeatMap = datasetService.create(new FloatType(), new long[]{dataset1.dimension(Axes.TIME), radialProfile.oCorrMap.keySet().size(), 3}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y, Axes.CHANNEL});
                    correlationAccessor = tempHeatMap.randomAccess();

                    ((LinearAxis)tempHeatMap.axis(1)).setScale(radialProfile.oCorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).max().getAsDouble()/radialProfile.oCorrMap.keySet().size());
                }

                double[] keySet = radialProfile.oCorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).toArray();

                for (int k = 0; k < keySet.length; k++) {
                    correlationAccessor.setPosition(new long[]{i,k,2});
                    correlationAccessor.get().setReal(radialProfile.oCorrMap.get(radialProfile.getBD(keySet[k])));
                    correlationAccessor.setPosition(new long[]{i,k,1});
                    correlationAccessor.get().setReal(radialProfile.sCorrMap.get(radialProfile.getBD(keySet[k])));
                    correlationAccessor.setPosition(new long[]{i,k,0});
                    correlationAccessor.get().setReal(radialProfile.gaussCurveMap.get(radialProfile.getBD(keySet[k])));
                }

                if(radialProfile.confidence > highestConfidence){
                    highestConfidence = radialProfile.confidence;
                    highestConMean = radialProfile.gaussFit[1];
                    highestConSD = radialProfile.gaussFit[2];
                    highestConFrame = i;
                    highestCCvalue = radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFit[1]));
                }

                LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<>();
                gaussianMap.put("Mean",  getSigDigits(radialProfile.gaussFit[1]));
                gaussianMap.put("SD",  getSigDigits(radialProfile.gaussFit[2]));
                gaussianMap.put("Gaussian height", getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFit[1]))));
                gaussianMap.put("Confidence",  getSigDigits(radialProfile.confidence));

                correlationTableList.add(gaussianMap);

                correlationTablesRowNames.add((calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "" + calibratedTime.get().calibratedValue(i) : "Frame " + i));
            }

            max = tempHeatMap.dimensionsAsLongArray();

            max[1] = Math.min(Math.round(tempHeatMap.axis(1).rawValue(highestConMean + 5 * highestConSD)), max[1]-1);

            RandomAccessibleInterval temp = ops.transform().crop(tempHeatMap, Intervals.createMinMax(0, 0, 0, max[0]-1, max[1], max[2]-1));

            timeCorrelationHeatMap = datasetService.create(temp);

            ((LinearAxis) timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1) : 1);
            timeCorrelationHeatMap.setAxis(tempHeatMap.axis(1).copy(), 1);

            timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit() : "frame"));
            timeCorrelationHeatMap.axis(0).setType(Axes.X);
            timeCorrelationHeatMap.axis(1).setUnit(getUnitType());
            timeCorrelationHeatMap.axis(1).setType(Axes.Y);
            //timeCorrelationHeatMap.axis(2).setType(Axes.CHANNEL);

            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Time of best CC (" + timeCorrelationHeatMap.axis(0).unit() + ")");
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Gaussian Height");
            rowHeaders.add("Confidence");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? getSigDigits(calibratedTime.get().calibratedValue(highestConFrame)) : highestConFrame);
            resultsList.add(getSigDigits(highestConMean));
            resultsList.add(getSigDigits(highestConSD));
            resultsList.add(getSigDigits(highestCCvalue));
            resultsList.add(getSigDigits(highestConfidence));

            results = Tables.wrap(resultsList, null, rowHeaders);

            if(showIntermediates){
                for (int i = 0; i < 4; i++) {
                    uiService.show(intermediates[i]);
                }
            }

            timeCorrelationHeatMap.setName("Heat map of correlation over time between " + dataset1.getName() + " and " + dataset2.getName());

            correlationTable = Tables.wrap(correlationTableList, correlationTablesRowNames);

            if(highestConfidence < 15){
                notes = (notes == null ? "" : notes) + "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.\n\n\n";
            }

            if(saveFolder != null) {
                if (!saveFolder.exists() || !saveFolder.canWrite()) {
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);
                    datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + "\\" + ContributionOf1.getName() + ".tif", config);
                    datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + "\\" + ContributionOf2.getName() + ".tif", config);

                    datasetIOService.save(timeCorrelationHeatMap, saveFolder.getAbsolutePath() + "\\" + timeCorrelationHeatMap.getName(), config);

                    ioService.save(Tables.wrap(correlationTableList, correlationTablesRowNames), saveFolder.getAbsolutePath() + "\\" + "Gaussian fits over time.csv");

                    String summary = (notes == null ? "" : notes) + "Highest confidence fit of a gaussian curve to the cross-correlation of: \n\""+
                            dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\"\n using the mask \n\"" + (maskAbsent? "No mask selected" : maskDataset.getName()) +
                            "\"\nwas found at " + (calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "time " + calibratedTime.get().calibratedValue(highestConFrame) + ", ": "") + "frame " + highestConFrame +
                            ":\n\nMean (" + getUnitType() +"): " + getSigDigits(highestConMean) +
                            "\nStandard deviation: " + getSigDigits(highestConSD) +
                            "\nGaussian height: " + getSigDigits(highestCCvalue) +
                            "\nConfidence: " + getSigDigits(highestConfidence) +
                            "\n\n\nThe 3-channel heat map shows the (by channel): \n 1. Gaussian curve for each frame.\n 2. Subtracted correlation for each frame.\n 3. Original correlation for each frame.\n\nFor more details, please see the website: \nhttps://imagej.github.io/Colocalization_by_Cross_Correlation";

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Summary.txt"), summary, (Charset) null);

                    if (showIntermediates) {
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + "\\" + intermediate.getName() + ".tif", config);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion
    }

    private String getUnitType(){
        return dataset1.axis(Axes.X).isPresent() ? dataset1.axis(Axes.X).get().unit(): "Unlabeled distance unit";
    }

    private double getVoxelVolume(){
        double volume = 1;
        for (int i = 0; i < scale.length; i++) {
            volume *= scale[i];
        }
        return volume;
    }

    //made this to quickly and easily test different extension methods for correlation
    private RandomAccessible extendImage(Img in){
        //return Views.extendMirrorSingle(in); //this is the default method, it causes major issues when there is a flat uniform background (even small numbers) over the whole image with no mask
        //return Views.extendValue(in, ops.stats().median(in).getRealFloat()); //this can cause issues similar to extendMirrorSingle, though slightly less often
        return Views.extendZero(in); //this method seems to be the best for cross-correlation. The original cross-correlation can look terrible with flat background or noise (looks like a pyramid), but this is subtracted out. This method also makes the most intuitive sense, as we don't want to correlate beyond the borders of the image.
    }

    private <T extends RealType> void colocalizationAnalysis(Img <? extends T> img1, Img <? extends T> img2, Img <? extends T> imgMask, RadialProfiler radialProfiler, final RandomAccessibleInterval <? extends RealType> contribution1, final RandomAccessibleInterval <? extends RealType> contribution2, RandomAccessibleInterval <? extends RealType> [] localIntermediates){
        statusService.showStatus(statusBase + "Applying masks");

        //Zero all the data outside the image mask, to prevent it from contributing to the cross-correlation result.
        LoopBuilder.setImages(img1, imgMask).multiThreaded().forEachPixel((a,b) -> {if((b.getRealFloat() == 0.0)) {a.setReal(b.getRealFloat());}});
        LoopBuilder.setImages(img2, imgMask).multiThreaded().forEachPixel((a,b) -> {if((b.getRealFloat() == 0.0)) {a.setReal(b.getRealFloat());}});

        statusService.showStatus(statusBase + "Initializing randomizer");
        CostesRandomizer imageRandomizer = new CostesRandomizer(img1, imgMask);

        double maskVolume = imageRandomizer.getMaskVoxelCount()*getVoxelVolume();

        statusService.showStatus(statusBase + "Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> oCorr = imgFactory.create(img1);
        Img<FloatType> rCorr = imgFactory.create(img1);

        final T zero = Util.getTypeFromInterval(img1);
        zero.setZero();

        OutOfBoundsFactory outOfBoundsFactory = new OutOfBoundsConstantValueFactory(zero);

        ops.filter().correlate(oCorr, img1, img2, img1.dimensionsAsLongArray(), outOfBoundsFactory);

        //normalize correlation product to mask volume
        LoopBuilder.setImages(oCorr).multiThreaded().forEachPixel((a) -> a.setReal(a.get()/maskVolume));


        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], oCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }


        /**Start creating average correlation of Costes Randomization data. Have to begin this outside the loop to seed
         * avgRandCorr with non-zero data. The zeroed data outside the mask is unaltered during the randomization process,
         * so that it does not contribute to the result. After we initialize it, we can continue to create an average
         * correlation with random data.
         *
         * While working with test data, I noticed that the number of randomizations is not crucial and often a single
         * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
         * data may require more randomization cycles.
         */

        Img <RealType> randomizedImage = imageRandomizer.getRandomizedImage();

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[1], randomizedImage).multiThreaded().forEachPixel((a,b) -> a.setReal(b.getRealFloat()));
        }


        //Need to change this, it weighs later cycles heavier than first ones. Change to add rCorr*1/cycles multiple times together
        statusService.showStatus(statusBase + "Cycle 1/" + cycles + " - Calculating randomized correlation");

        Img<FloatType> avgRandCorr = imgFactory.create(rCorr);

        for (int i = 0; i < cycles; ++i) {
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Randomizing Image");
            randomizedImage = imageRandomizer.getRandomizedImage();
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Calculating randomized correlation");
            ops.filter().correlate(rCorr, randomizedImage, img2, img1.dimensionsAsLongArray(), outOfBoundsFactory);
            LoopBuilder.setImages(rCorr).multiThreaded().forEachPixel((a) -> a.setReal(a.get()/maskVolume));
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Averaging randomized correlation");

            ImgMath.compute(ImgMath.add(rCorr, avgRandCorr)).into(avgRandCorr);
        }
        ImgMath.compute(ImgMath.div(avgRandCorr, cycles)).into(avgRandCorr);
        /*Subtract the random correlation from the original to generate a subtracted correlation map. This is
          what will be used to evaluate any spatial relations between the two channels.
         */

        statusService.showStatus(statusBase + "Subtracting randomized correlation");
        Img<FloatType> subtracted = oCorr.copy();
        ImgMath.compute(ImgMath.sub(oCorr, avgRandCorr)).into(subtracted);

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[2], subtracted).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }

        /* Plot the subtracted correlation in the same plot as the original data. Contribution from random elements
         * in the original data should lie close to zero relative to the original data. Real associations will be
         * less affected
         */
        statusService.showStatus(statusBase + "Calculating radial profile");
        try{radialProfiler.calculateProfiles(oCorr, subtracted);}
        catch (Exception e){
            DialogPrompt.Result result = uiService.showDialog("Failed to fit gaussian curve to data, suggesting no correlation between the images.\nSelect OK to show intermediate correlation images (if the option was selected). Select cancel to interrupt plugin and show full error message.", DialogPrompt.MessageType.ERROR_MESSAGE, DialogPrompt.OptionType.OK_CANCEL_OPTION);
            if (result == DialogPrompt.Result.CANCEL_OPTION){
                throw e;
            }
            else{
                for (RandomAccessibleInterval i:localIntermediates
                     ) {
                    uiService.show(i);
                }
                throw e;
            }
        }

        /*After getting the radial profile, need to fit a gaussian curve to the data, and draw the points to
         * the plot window.
         */

        /* Once we have the fit, we need to establish a confidence value in it. This is very important as the
         * gaussian fitter will always return a result, and even completely non-correlated images can occasionally
         * return low sigma values. The confidence is simply the area under the curve (range: mean +/- 3 sigma)
         * of the subtracted data curve divided by the area under the curve of the original data curve as a percentage. Particles with
         * strong spatial association are less affected by the subtraction of the randomized data and thus will
         * have higher confidence values.
         *
         * Decent values of confidence: ~>10 to >20?
         *
         * Note: It may be desirable to somehow scale the confidence based on the distance, as the confidence will
         * naturally decrease with increasing distance.
         */

        /* I finally figured it out! To get a representation of the signal from each image that contributed to the
         * cross-correlation after subtraction, I need to do a convolution between the subtracted correlation and
         * img(1?) , then multiply the result of that with the other image.
         *
         * Using rCorr for intermediate steps to avoid generating unnecessary Images
         */

        /*
        To get contributions, I can't just use subtracted, as it is equivalent to a lower valued version of the original
        correlation map. I have to modify subtracted by the Gaussian fit results
         */

        Img<FloatType> gaussModifiedCorr = subtracted.copy();
        ApplyGaussToCorr(subtracted, scale, radialProfiler.gaussCurveMap, gaussModifiedCorr, radialProfiler);
        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[3], gaussModifiedCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }
        statusService.showStatus(statusBase + "Determining channel contributions");

        //To get contribution of img1, convolve img2 with the gauss-modified correlation, then multiply with img1
        ops.filter().convolve(rCorr, img2, gaussModifiedCorr, img1.dimensionsAsLongArray(), outOfBoundsFactory);

        LoopBuilder.setImages(contribution1, ImgMath.compute(ImgMath.mul(rCorr, img1)).into(rCorr.copy())).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));

        //To get contribution of img2, correlate img1 with the gauss-modified correlation, then multiply with img2
        ops.filter().correlate(rCorr, img1, gaussModifiedCorr, img1.dimensionsAsLongArray(), outOfBoundsFactory);

        LoopBuilder.setImages(contribution2, ImgMath.compute(ImgMath.mul(rCorr, img2)).into(rCorr.copy())).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));

        //service.shutdown();
    }

    private double getSigDigits(double input){ return ((Math.round(input* sigDigits))/ sigDigits);}


    //troubleshooting method for showing images at key points
/*    private void showScaledImg(Img input, String title){
        uiService.show(title, input);
    }*/


    private <T extends RealType> void ApplyGaussToCorr(RandomAccessibleInterval <T> input, double[] scale, Map<BigDecimal, Double> gaussMap, RandomAccessibleInterval <T> output, RadialProfiler radialProfile){
        //get image dimensions and center
        int nDims = input.numDimensions();
        if(nDims != scale.length)
            return;
        long [] dims = new long[nDims];
        input.dimensions(dims);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dims[i])/2;
        }

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(input, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor <T> looper = Views.interval(input,chunk).localizingCursor();
                RandomAccess <T> outLooper = output.randomAccess();
                while(looper.hasNext()){
                    looper.fwd();
                    outLooper.setPosition(looper);
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i)-center[i])*scale[i],2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    outLooper.get().setReal(looper.get().getRealFloat()*gaussMap.get(radialProfile.getBD(Ldistance)));
                }
            });
        });
    }
}

