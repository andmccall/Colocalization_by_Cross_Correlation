import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.LinearAxis;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccess;
import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.math.ImgMath;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;
import net.imglib2.util.Intervals;
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
import org.scijava.util.ColorRGB;

import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, menuPath = "Analyze>Colocalization>Just Cross Correlation")
public class Just_Cross_Correlation implements Command{

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

    @Parameter(label = "Image 1: ", description = "This is the image which will be randomized during Costes randomization", persist = false)
    private Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    private Dataset dataset2;

    @Parameter(label = "Significant digits: ")
    private int significantDigits;

    @Parameter(label = "Output directory (leave blank for none):", description = "The directory to automatically save all generated output, including the intermediate images if the \"Show Intermediates\" box is checked", required = false, style="directory")
    private File saveFolder;

    @Parameter(type = ItemIO.OUTPUT)
    private Dataset timeCorrelationHeatMap;

    @Parameter(type = ItemIO.OUTPUT, label = "Table of correlations")
    private Table correlationTable;

    @Parameter(type = ItemIO.OUTPUT, label = "CC Results")
    private Table results;

    @Parameter(type = ItemIO.OUTPUT)
    private XYPlot plot;

    @Parameter(type = ItemIO.OUTPUT, label = "Summary of Results")
    private String notes;

    private double sigDigits;

    public SortedMap<BigDecimal, Double> CorrMap;

    private double [] scale;

    private String statusBase = "";

    public Just_Cross_Correlation() {
    }

    @Override
    public void run(){

        //region Error checking
        if(dataset1.numDimensions() != dataset2.numDimensions() || dataset1.getHeight() != dataset2.getHeight() || dataset1.getWidth() != dataset2.getWidth() || dataset1.getDepth() != dataset2.getDepth() || dataset1.getFrames() != dataset2.getFrames()){
            logService.error("All image dimensions (XYZ, and time) must match");
            return;
        }

        if(dataset1.getChannels() > 1 || dataset2.getChannels() > 1){
            logService.error("Multi-channel images are not supported, requires separate channels");
            return;
        }
        //endregion

        statusService.showStatus("Initializing plugin data");

        //region Plugin initialization (mostly creating datasets)
        RadialProfiler radialProfile;
        CorrMap = new TreeMap<>();
        sigDigits = Math.pow(10.0, significantDigits);

        SCIFIOConfig config = new SCIFIOConfig();
        config.writerSetFailIfOverwriting(false);

        // Cannot use duplicateBlank() for creating the upcoming images, as they need to be 32-bit Float images
        CalibratedAxis [] calibratedAxes = new CalibratedAxis[dataset1.numDimensions()];
        AxisType [] axisTypes = new AxisType[dataset1.numDimensions()];
        for (int i = 0; i < dataset1.numDimensions(); ++i) {
            calibratedAxes[i] = dataset1.axis(i);
            axisTypes[i] = dataset1.axis(i).type();
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
            try{colocalizationAnalysis(dataset1.duplicate(), dataset2.duplicate(),radialProfile);}
            catch (Exception e){
                e.printStackTrace();
                throw e;
            }


            SeriesStyle oCorrStyle = plotService.newSeriesStyle(ColorRGB.fromHTMLColor("blue"), LineStyle.SOLID, MarkerStyle.NONE);

            plot = plotService.newXYPlot();

            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(CorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).boxed().collect(Collectors.toList()), new ArrayList<>(CorrMap.values()));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");

            plot.xAxis().setLabel("Distance (" + getUnitType() + ")");

            plot.yAxis().setLabel("Cross correlation");

            //plot.xAxis().setManualRange(0, Math.max(radialProfile.gaussFitPamameters[1] + (5 * radialProfile.gaussFitPamameters[2]), 0.01));

            plot.setTitle("Correlation of images");

/*            if(radialProfile.confidence < 15){
                notes = (notes == null ? "" : notes) + "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.\n\n\n";
            }*/

/*            if(radialProfile.rSquared < 0.05){
                notes = (notes == null ? "" : notes) + "The R-squared value for the gaussian regression is low.\nThis usually indicates a low signal to noise ratio, and additional pre-processing steps may help, such as a";
            }*/

/*            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Gaussian Height");
            rowHeaders.add("Confidence");
            rowHeaders.add("R-squared");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(getSigDigits(radialProfile.gaussFitPamameters[1]));
            resultsList.add(getSigDigits(radialProfile.gaussFitPamameters[2]));
            resultsList.add(getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFitPamameters[1]))));
            resultsList.add(getSigDigits(radialProfile.confidence));
            resultsList.add(getSigDigits(radialProfile.rSquared));

            results = Tables.wrap(resultsList, "", rowHeaders);*/

            List<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            CorrMap.keySet().stream().forEachOrdered((d) -> {
                LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
                row.put("Distance (" + getUnitType() +")", (getSigDigits(d.doubleValue())));
                row.put("Cross-correlation", getSigDigits(CorrMap.get(d)));
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

                    File plotout = new File(saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".png");
                    XYPlotConverter converter = new XYPlotConverter();
                    ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth()*2, plot.getPreferredHeight()*2);

                    //ioService.save(results,saveFolder.getAbsolutePath() + "\\" + "CC Results.csv" );
                    ioService.save(correlationTable, saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".csv");

/*                    String summary = (notes == null ? "" : notes) + "Fit a gaussian curve to the cross-correlation of: \n\""
                            + dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\":\n\nMean (" + getUnitType() +"): " + getSigDigits(radialProfile.gaussFitPamameters[1]) +
                            "\nStandard deviation: " + getSigDigits(radialProfile.gaussFitPamameters[2]) +
                            "\nGaussian height:" + getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFitPamameters[1])));

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Summary.txt"), summary, (Charset) null);*/

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

/*            ArrayList<HashMap<String,Double>> correlationTableList = new ArrayList<>();
            ArrayList<String> correlationTablesRowNames = new ArrayList<>();

            double highestConfidence = 0;
            double highestConMean = 0;
            double highestConSD = 0;
            long highestConFrame = 0;
            double highestCCvalue = 0;
            double highestRsquared = 0;*/


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

                try {
                    radialProfile = new RadialProfiler(temp1, scale);
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }

                try{colocalizationAnalysis(datasetService.create(temp1), datasetService.create(temp2), radialProfile);}
                catch (Exception e){
                    e.printStackTrace();
                    throw e;
                }

                if(tempHeatMap == null) {
                    tempHeatMap = datasetService.create(new FloatType(), new long[]{dataset1.dimension(Axes.TIME), CorrMap.keySet().size()}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y});
                    correlationAccessor = tempHeatMap.randomAccess();

                    ((LinearAxis)tempHeatMap.axis(1)).setScale(CorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).max().getAsDouble()/CorrMap.keySet().size());
                }

                double[] keySet = CorrMap.keySet().stream().mapToDouble(BigDecimal::doubleValue).toArray();

                for (int k = 0; k < keySet.length; k++) {
                    correlationAccessor.setPosition(new long[]{i,k});
                    correlationAccessor.get().setReal(CorrMap.get(radialProfile.getBD(keySet[k])));
                }

/*                if(radialProfile.confidence > highestConfidence){
                    highestConfidence = radialProfile.confidence;
                    highestConMean = radialProfile.gaussFitPamameters[1];
                    highestConSD = radialProfile.gaussFitPamameters[2];
                    highestConFrame = i;
                    highestCCvalue = radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFitPamameters[1]));
                    highestRsquared = radialProfile.rSquared;
                }

                LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<>();
                gaussianMap.put("Mean",  getSigDigits(radialProfile.gaussFitPamameters[1]));
                gaussianMap.put("SD",  getSigDigits(radialProfile.gaussFitPamameters[2]));
                gaussianMap.put("Gaussian height", getSigDigits(radialProfile.gaussCurveMap.get(radialProfile.getBD(radialProfile.gaussFitPamameters[1]))));
                gaussianMap.put("Confidence",  getSigDigits(radialProfile.confidence));
                gaussianMap.put("R-squared", getSigDigits(radialProfile.rSquared));

                correlationTableList.add(gaussianMap);



                correlationTablesRowNames.add((calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "" + calibratedTime.get().calibratedValue(i) : "Frame " + i));
 */
            }


            timeCorrelationHeatMap = tempHeatMap.copy();

            ((LinearAxis) timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1) : 1);
            timeCorrelationHeatMap.setAxis(tempHeatMap.axis(1).copy(), 1);

            timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit() : "frame"));
            timeCorrelationHeatMap.axis(0).setType(Axes.X);
            timeCorrelationHeatMap.axis(1).setUnit(getUnitType());
            timeCorrelationHeatMap.axis(1).setType(Axes.Y);
            //timeCorrelationHeatMap.axis(2).setType(Axes.CHANNEL);

/*            List<String> rowHeaders = new ArrayList<>();
            rowHeaders.add("Time of best CC (" + timeCorrelationHeatMap.axis(0).unit() + ")");
            rowHeaders.add("Mean (" + getUnitType() + ")");
            rowHeaders.add("StDev (" + getUnitType() + ")");
            rowHeaders.add("Gaussian Height");
            rowHeaders.add("Confidence");
            rowHeaders.add("R-squared");

            List<Double> resultsList = new ArrayList<>();
            resultsList.add(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? getSigDigits(calibratedTime.get().calibratedValue(highestConFrame)) : highestConFrame);
            resultsList.add(getSigDigits(highestConMean));
            resultsList.add(getSigDigits(highestConSD));
            resultsList.add(getSigDigits(highestCCvalue));
            resultsList.add(getSigDigits(highestConfidence));
            resultsList.add(getSigDigits(highestRsquared));

            results = Tables.wrap(resultsList, null, rowHeaders);*/

            timeCorrelationHeatMap.setName("Heat map of correlation over time between " + dataset1.getName() + " and " + dataset2.getName());

            //correlationTable = Tables.wrap(correlationTableList, correlationTablesRowNames);

/*            if(highestConfidence < 15){
                notes = (notes == null ? "" : notes) + "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.\n\n\n";
            }*/

            if(saveFolder != null) {
                if (!saveFolder.exists() || !saveFolder.canWrite()) {
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);

                    datasetIOService.save(timeCorrelationHeatMap, saveFolder.getAbsolutePath() + "\\" + timeCorrelationHeatMap.getName(), config);
/*
                    ioService.save(Tables.wrap(correlationTableList, correlationTablesRowNames), saveFolder.getAbsolutePath() + "\\" + "Gaussian fits over time.csv");

                    String summary = (notes == null ? "" : notes) + "Highest confidence fit of a gaussian curve to the cross-correlation of: \n\""+
                            dataset1.getName() +
                            "\"\n with \n\"" +
                            dataset2.getName() +
                            "\"\nwas found at " + (calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "time " + calibratedTime.get().calibratedValue(highestConFrame) + ", ": "") + "frame " + highestConFrame +
                            ":\n\nMean (" + getUnitType() +"): " + getSigDigits(highestConMean) +
                            "\nStandard deviation: " + getSigDigits(highestConSD) +
                            "\nGaussian height: " + getSigDigits(highestCCvalue) +
                            "\nConfidence: " + getSigDigits(highestConfidence) +
                            "\nR-squared: " + getSigDigits(highestRsquared) +
                            "\n\n\nThe 3-channel heat map shows the (by channel): \n 1. Gaussian curve for each frame.\n 2. Subtracted correlation for each frame.\n 3. Original correlation for each frame.\n\nFor more details, please see the website: \nhttps://imagej.github.io/Colocalization_by_Cross_Correlation";

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Summary.txt"), summary, (Charset) null);
*/

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

    //made this to quickly and easily test different extension methods for correlation
    private RandomAccessible extendImage(Img in){
        //return Views.extendMirrorSingle(in); //this is the default method, it causes major issues when there is a flat uniform background (even small numbers) over the whole image with no mask
        //return Views.extendValue(in, ops.stats().median(in).getRealDouble()); //this can cause issues similar to extendMirrorSingle, though slightly less often
        return Views.extendZero(in); //this method seems to be the best for cross-correlation. The original cross-correlation can look terrible with flat background or noise (looks like a pyramid), but this is subtracted out. This method also makes the most intuitive sense, as we don't want to correlate beyond the borders of the image.
    }

    private <T extends RealType> void colocalizationAnalysis(Img <? extends T> img1, Img <? extends T> img2, RadialProfiler radialProfiler){
               statusService.showStatus(statusBase + "Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> oCorr = imgFactory.create(img1);

        //OutOfBoundsFactory zeroBounds = new OutOfBoundsConstantValueFactory<>(0.0);
        //ops.filter().correlate(oCorr, img1, img2, img1.dimensionsAsLongArray(), zeroBounds, zeroBounds);



        ExecutorService service = Executors.newCachedThreadPool();

        FFTConvolution conj = new FFTConvolution(extendImage(img1), img1, Views.extendZero(img2), img2, img1.factory().imgFactory( new ComplexFloatType() ),  service);
        conj.setComputeComplexConjugate(true);
        conj.setOutput(oCorr);
        conj.convolve();

        //normalize correlation product to mask volume
        //LoopBuilder.setImages(oCorr).multiThreaded().forEachPixel((a) -> a.setReal(a.get()/maskVolume));

        statusService.showStatus(statusBase + "Calculating radial profile");
        radialProfiler.calculateSingleProfile(oCorr, CorrMap);

        service.shutdown();
    }

    private double getSigDigits(double input){ return ((Math.round(input* sigDigits))/ sigDigits);}


    //troubleshooting method for showing images at key points
/*    private void showScaledImg(Img input, String title){
        uiService.show(title, input);
    }*/

}

