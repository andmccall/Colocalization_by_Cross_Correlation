import io.scif.config.SCIFIOConfig;

import io.scif.services.DatasetIOService;
import net.imagej.*;

import net.imagej.axis.*;
import net.imagej.display.ColorTables;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;

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
import org.scijava.ui.UIService;
import org.scijava.ui.swing.viewer.plot.jfreechart.XYPlotConverter;
import org.scijava.util.*;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.nio.charset.Charset;
import java.util.*;
import static java.util.stream.Collectors.toList;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, headless = true, menuPath = "Analyze>Colocalization>Colocalization by Cross Correlation")
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
    private PlotService plotService;

    @Parameter
    private IOService ioService;

    @Parameter
    private DatasetIOService datasetIOService;

    @Parameter
    private OpService ops;

    @Parameter(label = "Image 1: ", description = "This is the image which will be randomized during pixel randomization", persist = false)
    private Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    private Dataset dataset2;

    @Parameter(label = "No mask (not recommended)?", description = "When checked, performs pixel randomization over the entire image, regardless of what image is selected below.", callback = "maskCallback")
    private boolean maskAbsent;

    @Parameter(label = "Mask: ", description = "The mask over which pixels of image 1 will be randomized. This is important, more details at: imagej.github.io/Colocalization_by_Cross_Correlation", required = false, persist = false)
    private Dataset maskDataset;

    @Parameter(label = "Significant digits: ")
    private int significantDigits;

    @Parameter(label = "Generate contribution images?", description = "Generates images that highlight the signal from Image 1 and Image 2 that contributed to the result. Uncheck to use less memory.")
    private boolean generateContributionImages;

    @Parameter(label = "Show intermediate images? ", description = "Shows images of numerous steps throughout the algorithm. Uncheck to use less memory. More details at: imagej.github.io/Colocalization_by_Cross_Correlation")
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

    private Dataset [] intermediates;

    private String [] intermediateNames = {"Original CC result", "Subtracted CC result", "Gaussian-modified CC result"};

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
        String version = "2.0.0";
        String summary = "Results generated using CCC version " + version + "\n"
                + "Plugin website: https://imagej.net/plugins/colocalization-by-cross-correlation\n\n";

        if(saveFolder != null) saveFolder.mkdirs();

        RadialProfiler radialProfile = null;

        if(maskAbsent){
            maskDataset = dataset1.duplicateBlank();
            LoopBuilder.setImages(maskDataset).multiThreaded().forEachPixel(SetOne::setOne);
        }

        SCIFIOConfig config = new SCIFIOConfig();
        config.writerSetFailIfOverwriting(false);

        {
            // Cannot use duplicateBlank() for creating the upcoming images, as they need to be 32-bit Float images
            CalibratedAxis[] calibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
            AxisType[] axisTypes = new AxisType[dataset1.numDimensions()];
            for (int i = 0; i < dataset1.numDimensions(); ++i) {
                calibratedAxes[i] = dataset1.axis(i);
                axisTypes[i] = dataset1.axis(i).type();
            }

            if (generateContributionImages) {
                ContributionOf1 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset1.getName(), axisTypes);
                ContributionOf1.setAxes(calibratedAxes);
                ContributionOf2 = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), "Contribution of " + dataset2.getName(), axisTypes);
                ContributionOf2.setAxes(calibratedAxes);
            }

            if (showIntermediates) {
                intermediates = new Dataset[3];

                for (int i = 0; i < 3; i++) {
                    intermediates[i] = datasetService.create(new FloatType(), dataset1.dimensionsAsLongArray(), intermediateNames[i], axisTypes);
                    intermediates[i].setAxes(calibratedAxes);
                }
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

            if(maskAbsent){
                try {
                    colocalizationAnalysis(dataset1, dataset2, maskDataset, radialProfile, ContributionOf1, ContributionOf2, intermediates);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            //Have to duplicate the datasets to be able to apply the mask without changing the original input images
            else {
                try {
                    colocalizationAnalysis(dataset1.duplicate(), dataset2.duplicate(), maskDataset, radialProfile, ContributionOf1, ContributionOf2, intermediates);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }

            if(showIntermediates){
                for (Dataset intermediate : intermediates) {
                    uiService.show(intermediate);
                }
            }

            SeriesStyle oCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("blue"), LineStyle.NONE, MarkerStyle.CIRCLE);
            SeriesStyle sCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("green"), LineStyle.NONE, MarkerStyle.CIRCLE);
            SeriesStyle gaussStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("#ff00ff"), LineStyle.NONE, MarkerStyle.FILLEDCIRCLE);

            plot = plotService.newXYPlot();

            XYSeries gaussData = plot.addXYSeries();
            gaussData.setValues(new ArrayList<>(radialProfile.sCorrMap.keySet()), radialProfile.sCorrMap.keySet().stream().map(radialProfile.gaussian::value).collect(toList()));
            gaussData.setStyle(gaussStyle);
            gaussData.setLabel("Gaussian Fit");

            XYSeries sCorrPlotData = plot.addXYSeries();
            sCorrPlotData.setValues(new ArrayList<>(radialProfile.sCorrMap.keySet()), new ArrayList<>(radialProfile.sCorrMap.values()));
            sCorrPlotData.setStyle(sCorrStyle);
            sCorrPlotData.setLabel("Subtracted CC");

            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(new ArrayList<>(radialProfile.oCorrMap.keySet()), new ArrayList<>(radialProfile.oCorrMap.values()));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");

            plot.xAxis().setLabel("Distance (" + getUnitType() + ")");

            plot.yAxis().setLabel("Cross correlation");

            plot.xAxis().setManualRange(0, Math.max(radialProfile.gaussFitParameters[1] + (5 * radialProfile.gaussFitParameters[2]), 0.01));

            plot.setTitle("Correlation of images");

            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Confidence");
            rowHeaders.add("R-squared");
            rowHeaders.add("Gaussian Height");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(getSigDigits(radialProfile.gaussFitParameters[1]));
            resultsList.add(getSigDigits(radialProfile.gaussFitParameters[2]));
            resultsList.add(getSigDigits(radialProfile.confidence));
            resultsList.add(getSigDigits(radialProfile.rSquared));
            resultsList.add(getSigDigits(radialProfile.gaussian.value(radialProfile.gaussFitParameters[1])));

            results = Tables.wrap(resultsList, "", rowHeaders);

            List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            RadialProfiler finalRadialProfile = radialProfile;
            radialProfile.oCorrMap.keySet().stream().forEachOrdered((d) -> {
                LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
                row.put("Distance (" + getUnitType() +")", (getSigDigits(d)));
                row.put("Original CC", getSigDigits(finalRadialProfile.oCorrMap.get(d)));
                row.put("Subtracted CC", getSigDigits(finalRadialProfile.sCorrMap.get(d)));
                row.put("Gaussian fit", getSigDigits(finalRadialProfile.gaussian.value(d)));
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

                    File plotout = new File(saveFolder.getAbsolutePath() + File.separator + plot.getTitle() + ".png");
                    XYPlotConverter converter = new XYPlotConverter();
                    ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth()*2, plot.getPreferredHeight()*2);

                    ioService.save(results,saveFolder.getAbsolutePath() + File.separator + "CC Results.csv" );
                    ioService.save(correlationTable, saveFolder.getAbsolutePath() + File.separator + plot.getTitle() + ".csv");

                    summary = summary + "Fit a gaussian curve to the cross-correlation of: \n\""
                            + dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\"\n using the mask \n\"" +
                            (maskAbsent? "No mask selected" : maskDataset.getName()) +
                            "\":\n\nMean (" + getUnitType() +"): " + getSigDigits(radialProfile.gaussFitParameters[1]) +
                            "\nStandard deviation: " + getSigDigits(radialProfile.gaussFitParameters[2]) +
                            "\n\nConfidence: " + getSigDigits(radialProfile.confidence) +
                            "\nR-squared: " + getSigDigits(radialProfile.rSquared) +
                            "\n\nGaussian height (generally unused):" + getSigDigits(radialProfile.gaussian.value(radialProfile.gaussFitParameters[1]));

                            if(radialProfile.confidence < 0.1){
                                summary = summary + "\n\nThe confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nSee the website for more details.";
                            }

                            if(radialProfile.rSquared < 0.05){
                                summary = summary + "\n\nThe R-squared value for the gaussian regression is very low.\nThis can indicate a low signal to noise ratio, or that no spatial correlation exists and the curve was fit to image noise.\nSee website for more details.";
                            }

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + File.separator + "Summary.txt"), summary, (Charset) null);

                    if(generateContributionImages) {
                        datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + File.separator + ContributionOf1.getName() + ".tif", config);
                        datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + File.separator + ContributionOf2.getName() + ".tif", config);
                    }

                    if(showIntermediates){
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + File.separator + intermediate.getName() + ".tif", config);
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
            double highestRsquared = 0;


            RandomAccessibleInterval [] intermediatesViewsPasser = null;
            if(showIntermediates){
                intermediatesViewsPasser = new RandomAccessibleInterval[intermediates.length];
            }

            //Duplicate datasets to not modify originals when applying masks later
            Dataset dataset1copy = dataset1.duplicate();
            Dataset dataset2copy = dataset2.duplicate();

            //Dataset tempHeatMap = null;

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
                    for (int m = 0; m < intermediates.length; m++) {
                        intermediatesViewsPasser[m] = Views.dropSingletonDimensions(Views.interval(intermediates[m], min, max));
                    }
                }

                try {
                    radialProfile = new RadialProfiler(temp1, scale);
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }


                try {
                    colocalizationAnalysis(datasetService.create(temp1), datasetService.create(temp2), datasetService.create(masktemp), radialProfile, ContributionOf1 == null ? null : Views.dropSingletonDimensions(Views.interval(ContributionOf1, min, max)), ContributionOf2 == null ? null : Views.dropSingletonDimensions(Views.interval(ContributionOf2, min, max)), intermediatesViewsPasser);
                } catch (Exception e) {
                    e.printStackTrace();
                    throw e;
                }


                if(timeCorrelationHeatMap == null) {
                    timeCorrelationHeatMap = datasetService.create(new FloatType(), new long[]{dataset1.dimension(Axes.TIME), radialProfile.oCorrMap.keySet().size(), 3}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y, Axes.CHANNEL});

                    ((LinearAxis) timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1) : 1);
                    timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit() : "frame"));
                    timeCorrelationHeatMap.axis(0).setType(Axes.X);
                    timeCorrelationHeatMap.axis(2).setType(Axes.CHANNEL);
                    timeCorrelationHeatMap.initializeColorTables(3);
                    timeCorrelationHeatMap.setColorTable(ColorTables.MAGENTA,0);

                    ((LinearAxis) timeCorrelationHeatMap.axis(1)).setScale((radialProfile.oCorrMap.lastKey()-radialProfile.oCorrMap.firstKey())/radialProfile.oCorrMap.keySet().size());
                    ((LinearAxis) timeCorrelationHeatMap.axis(1)).setOrigin(radialProfile.oCorrMap.firstKey());
                    //EnumeratedAxis seems to be broken and doesn't show the proper values in ImageJ
                    //timeCorrelationHeatMap.setAxis(new EnumeratedAxis(Axes.Y, getUnitType(), radialProfile.oCorrMap.keySet().stream().mapToDouble(Double::doubleValue).toArray()), 1);
                    //uiService.showDialog("First key: " + radialProfile.oCorrMap.firstKey() + "," + timeCorrelationHeatMap.axis(1).rawValue(radialProfile.oCorrMap.firstKey()));

                    correlationAccessor = timeCorrelationHeatMap.randomAccess();
                }

                double[] keySet = radialProfile.oCorrMap.keySet().stream().mapToDouble(Double::doubleValue).toArray();

                for (int k = 0; k < keySet.length; k++) {
                    correlationAccessor.setPosition(new long[]{i,k,2});
                    correlationAccessor.get().setReal(radialProfile.oCorrMap.get(keySet[k]));
                    correlationAccessor.setPosition(new long[]{i,k,1});
                    correlationAccessor.get().setReal(radialProfile.sCorrMap.get(keySet[k]));
                    correlationAccessor.setPosition(new long[]{i,k,0});
                    correlationAccessor.get().setReal(radialProfile.gaussian.value(keySet[k]));
                }

                if(radialProfile.confidence > highestConfidence){
                    highestConfidence = radialProfile.confidence;
                    highestConMean = radialProfile.gaussFitParameters[1];
                    highestConSD = radialProfile.gaussFitParameters[2];
                    highestConFrame = i;
                    highestCCvalue = radialProfile.gaussian.value(radialProfile.gaussFitParameters[1]);
                    highestRsquared = radialProfile.rSquared;
                }

                LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<>();
                gaussianMap.put("Mean",  getSigDigits(radialProfile.gaussFitParameters[1]));
                gaussianMap.put("SD",  getSigDigits(radialProfile.gaussFitParameters[2]));
                gaussianMap.put("Confidence",  getSigDigits(radialProfile.confidence));
                gaussianMap.put("R-squared", getSigDigits(radialProfile.rSquared));
                gaussianMap.put("Gaussian height", getSigDigits(radialProfile.gaussian.value(radialProfile.gaussFitParameters[1])));

                correlationTableList.add(gaussianMap);

                correlationTablesRowNames.add((calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "" + getSigDigits(calibratedTime.get().calibratedValue(i)) + " " + calibratedTime.get().unit() : "Frame " + i));
            }

            /*max = tempHeatMap.dimensionsAsLongArray();

            max[1] = Math.min(Math.round(tempHeatMap.axis(1).rawValue(highestConMean + 5 * highestConSD)), max[1]-1);

            RandomAccessibleInterval temp = ops.transform().crop(tempHeatMap, Intervals.createMinMax(0, 0, 0, max[0]-1, max[1], max[2]-1));*/

            //timeCorrelationHeatMap = timeCorrelationHeatMap.copy();

            //timeCorrelationHeatMap.setAxis(timeCorrelationHeatMap.axis(1).copy(), 1);

            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Time of best CC (" + timeCorrelationHeatMap.axis(0).unit() + ")");
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Confidence");
            rowHeaders.add("R-squared");
            rowHeaders.add("Gaussian Height");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? getSigDigits(calibratedTime.get().calibratedValue(highestConFrame)) : highestConFrame);
            resultsList.add(getSigDigits(highestConMean));
            resultsList.add(getSigDigits(highestConSD));
            resultsList.add(getSigDigits(highestConfidence));
            resultsList.add(getSigDigits(highestRsquared));
            resultsList.add(getSigDigits(highestCCvalue));

            results = Tables.wrap(resultsList, null, rowHeaders);

            if(showIntermediates){
                for (Dataset intermediate : intermediates) {
                    uiService.show(intermediate);
                }
            }

            timeCorrelationHeatMap.setName("Heat map of correlation over time between " + dataset1.getName() + " and " + dataset2.getName());

            correlationTable = Tables.wrap(correlationTableList, correlationTablesRowNames);


            if(saveFolder != null) {
                if (!saveFolder.exists() || !saveFolder.canWrite()) {
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);
                    datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + File.separator + ContributionOf1.getName() + ".tif", config);
                    datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + File.separator + ContributionOf2.getName() + ".tif", config);

                    datasetIOService.save(timeCorrelationHeatMap, saveFolder.getAbsolutePath() + File.separator + timeCorrelationHeatMap.getName(), config);

                    ioService.save(Tables.wrap(correlationTableList, correlationTablesRowNames), saveFolder.getAbsolutePath() + File.separator + "Gaussian fits over time.csv");

                    summary = summary + "Highest confidence fit of a gaussian curve to the cross-correlation of: \n\""+
                            dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\"\n using the mask \n\"" + (maskAbsent? "No mask selected" : maskDataset.getName()) +
                            "\"\nwas found at " + (calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "time " + calibratedTime.get().calibratedValue(highestConFrame) + ", ": "") + "frame " + highestConFrame +
                            ":\n\nMean (" + getUnitType() +"): " + getSigDigits(highestConMean) +
                            "\nStandard deviation: " + getSigDigits(highestConSD) +
                            "\n\nConfidence: " + getSigDigits(highestConfidence) +
                            "\nR-squared: " + getSigDigits(highestRsquared) +
                            "\n\nGaussian height (generally unused): " + getSigDigits(highestCCvalue) +
                            "\n\n\nThe 3-channel heat map shows the (by channel): \n 1. Gaussian curve for each frame.\n 2. Subtracted correlation for each frame.\n 3. Original correlation for each frame.\n\nFor more details, please see the website.";

                    if(highestConfidence < 0.1){
                        summary = summary + "\n\nThe confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nSee the website for more details.";
                    }

                    if(highestRsquared < 0.05){
                        summary = summary + "\n\nThe R-squared value for the gaussian regression is very low.\nThis can indicate a low signal to noise ratio, or that no spatial correlation exists and the curve was fit to image noise.\nSee website for more details.";
                    }

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + File.separator + "Summary.txt"), summary, (Charset) null);

                    if (showIntermediates) {
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + File.separator + intermediate.getName() + ".tif", config);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion
    }

    private String getUnitType(){ return dataset1.axis(Axes.X).isPresent() ? dataset1.axis(Axes.X).get().unit(): "Unlabeled distance unit"; }

    private double getSigDigits(double input){
        BigDecimal bd = new BigDecimal(input);
        bd = bd.round(new MathContext(significantDigits));
        return bd.doubleValue();
        //return ((Math.round(input* sigDigits))/ sigDigits);
    }

    private <T extends RealType> void colocalizationAnalysis(Img <T> img1, Img <T> img2, Img <T> imgMask, RadialProfiler radialProfiler, final RandomAccessibleInterval <T> contribution1, final RandomAccessibleInterval <T> contribution2, RandomAccessibleInterval <T> [] localIntermediates){
        Img<FloatType> oCorr = ops.create().img(img1, new FloatType());
        Img<FloatType> subtracted = ops.create().img(img1, new FloatType());
        Img<FloatType> gaussModifiedCorr;

        //Zero all the data outside the image mask, to prevent it from contributing to the cross-correlation result.
        if(!maskAbsent) {
            statusService.showStatus(statusBase + "Applying masks");
            LoopBuilder.setImages(img1, imgMask).multiThreaded().forEachPixel((a, b) -> {
                if ((b.getRealDouble() == 0.0)) {
                    a.setReal(b.getRealDouble());
                }
            });
            LoopBuilder.setImages(img2, imgMask).multiThreaded().forEachPixel((a, b) -> {
                if ((b.getRealDouble() == 0.0)) {
                    a.setReal(b.getRealDouble());
                }
            });
        }

        statusService.showStatus(statusBase + "Generating averaged mask");

        CCfunctions ccFunctions = new CCfunctions(img1, img2, imgMask, scale);

        statusService.showStatus(statusBase + "Calculating original correlation");

        ccFunctions.calculateCC(oCorr);

        statusService.showStatus(statusBase + "Generating subtracted correlation");

        ccFunctions.generateSubtractedCCImage(img1, img2, imgMask, oCorr, subtracted);

        statusService.showStatus(statusBase + "Calculating radial profile");
        radialProfiler.calculateProfiles(oCorr, subtracted);

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], oCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }
        oCorr = null;

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[1], subtracted).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }

        if(!generateContributionImages){
            subtracted = null;
        }

        statusService.showStatus(statusBase + "Fitting gaussian to data");
        try{radialProfiler.fitGaussianCurve();}
        catch (NullPointerException e){
            generateContributionImages = false;
            logService.warn("Failed to fit gaussian curve to cross correlation of " + dataset1.getName() + " and " + dataset2.getName() + ", suggesting no correlation between the images.\nAcquired data and intermediate correlation images (if the option was selected) will still be shown. Statistical measures will be set to error values (-1).");
        }

        if(generateContributionImages) {
            statusService.showStatus(statusBase + "Determining channel contributions");
            //gaussModifiedCorr = imgFactory.create(img1);
            gaussModifiedCorr = ops.create().img(img1, new FloatType());

            ccFunctions.generateGaussianModifiedCCImage(subtracted, gaussModifiedCorr, radialProfiler);

            if(showIntermediates){
                LoopBuilder.setImages(localIntermediates[2], gaussModifiedCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
            }

            ccFunctions.calculateContributionImages(img1, img2, gaussModifiedCorr, contribution1, contribution2);
        }
    }
}

