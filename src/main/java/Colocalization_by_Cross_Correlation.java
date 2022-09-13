

import io.scif.config.SCIFIOConfig;
import io.scif.services.DatasetIOService;
import net.imagej.*;

import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imagej.axis.CalibratedAxis;
import net.imagej.axis.LinearAxis;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.RandomAccess;
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
import org.scijava.table.Tables;
import org.scijava.ui.DialogPrompt;
import org.scijava.ui.UIService;
import org.scijava.ui.swing.viewer.plot.jfreechart.XYPlotConverter;
import org.scijava.util.*;

import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


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

    @Parameter(type = ItemIO.OUTPUT)
    private XYPlot plot;

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

        // Cannot use duplicateBlank() for creating the upcoming images, as they need to be 32-bit
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

            XYSeries oCorrPlotData = plot.addXYSeries();
            oCorrPlotData.setValues(Arrays.asList(radialProfile.Xvalues), Arrays.asList(radialProfile.Yvalues[0]));
            oCorrPlotData.setStyle(oCorrStyle);
            oCorrPlotData.setLabel("Original CC");

            XYSeries sCorrPlotData = plot.addXYSeries();
            sCorrPlotData.setValues(Arrays.asList(radialProfile.Xvalues), Arrays.asList(radialProfile.Yvalues[1]));
            sCorrPlotData.setStyle(sCorrStyle);
            sCorrPlotData.setLabel("Subtracted CC");

            XYSeries gaussData = plot.addXYSeries();
            gaussData.setValues(Arrays.asList(radialProfile.Xvalues), Arrays.asList(radialProfile.Yvalues[2]));
            gaussData.setStyle(gaussStyle);
            gaussData.setLabel("Gaussian Fit");

            plot.xAxis().setManualRange(0, radialProfile.gaussFit[1] + (5* radialProfile.gaussFit[2]));

            plot.setTitle("Correlation of images");

            float mean = (float) getSigDigits(radialProfile.gaussFit[1]);

            String output = "Fit a gaussian curve to the cross-correlation of: \n\""
                    + dataset1.getName() +
                    "\"\n with \n\"" +
                    dataset2.getName() +
                    "\"\n using the mask \n\"" +
                    (maskAbsent? "No mask selected" : maskDataset.getName()) +
                    "\":\n\nMean: " + mean +
                    "\nStandard deviation (sigma): " + getSigDigits(radialProfile.gaussFit[2]) +
                    "\nGaussian height:" + getSigDigits(radialProfile.Yvalues[2][Math.round(mean)]) +
                    "\nConfidence: " + getSigDigits(radialProfile.confidence);

            uiService.show("Gauss Fit", output);

            if(radialProfile.confidence < 15){
                uiService.show("Low confidence", "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.");
            }

            if(saveFolder != null){
                if(!saveFolder.exists() || !saveFolder.canWrite()){
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);
                    datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + "\\" + ContributionOf1.getName(), config);
                    datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + "\\" + ContributionOf2.getName(), config);

                    File plotout = new File(saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".png");
                    XYPlotConverter converter = new XYPlotConverter();

                   ChartUtils.saveChartAsPNG(plotout, converter.convert(plot, JFreeChart.class), plot.getPreferredWidth()*2, plot.getPreferredHeight()*2);

                   List listOfPlotPoints = new ArrayList<HashMap>();

                    for (int i = 0; i < radialProfile.Xvalues.length; i++) {
                        LinkedHashMap<String, Double> row = new LinkedHashMap<String, Double>();
                        row.put("X", radialProfile.Xvalues[i]);
                        row.put("Original CC", radialProfile.Yvalues[0][i]);
                        row.put("Subtracted CC", radialProfile.Yvalues[1][i]);
                        row.put("Gaussian fit", radialProfile.Yvalues[2][i]);
                        listOfPlotPoints.add(row);
                    }

                    ioService.save(Tables.wrap(listOfPlotPoints, null), saveFolder.getAbsolutePath() + "\\" + plot.getTitle() + ".csv");

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Results Summary.txt"), output, (Charset) null);

                    if(showIntermediates){
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + "\\" + intermediate.getName(), config);
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

            max[timeAxis] = 0;
            timeCorrelationHeatMap = datasetService.create(new FloatType(), new long []{dataset1.dimension(Axes.TIME), RadialProfiler.getNumberOfBins(Views.dropSingletonDimensions(Views.interval(dataset1, min, max)), scale), 3}, "Correlation over time of " + dataset1.getName() + " and " + dataset2.getName(), new AxisType[]{Axes.X, Axes.Y, Axes.CHANNEL} );
            RandomAccess<RealType<?>> correlationAccessor = timeCorrelationHeatMap.randomAccess();

            ((LinearAxis)timeCorrelationHeatMap.axis(0)).setScale(calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? calibratedTime.get().calibratedValue(1): 1);
            ((LinearAxis)timeCorrelationHeatMap.axis(1)).setScale(RadialProfiler.getBinSize(Views.dropSingletonDimensions(Views.interval(dataset1, min, max)), scale));
            timeCorrelationHeatMap.axis(0).setUnit((dataset1.axis(Axes.TIME).isPresent() ? dataset1.axis(Axes.TIME).get().unit(): "Unlabeled time unit"));
            timeCorrelationHeatMap.axis(1).setUnit((dataset1.axis(Axes.X).isPresent() ? dataset1.axis(Axes.X).get().unit(): "Unlabeled distance unit"));

            List listOfGaussianMaps = new ArrayList();
            List rowNames = new ArrayList();

            double highestConfidence = 0;
            double highestConMean = 0;
            double highestConSD = 0;
            long highestConFrame = 0;
            double highestCCvalue = 0;

            RandomAccessibleInterval <? extends RealType> [] intermediatesViewsPasser = null;
            if(showIntermediates){
                intermediatesViewsPasser = new RandomAccessibleInterval[4];
            }

            for (long i = 0; i < dataset1.getFrames(); i++) {
                statusService.showProgress((int)i, (int)dataset1.getFrames());
                //StatusService message update is performed in colocalizationAnalysis method
                statusBase = "Frame " + (i+1) + " - ";
                min[timeAxis] = i;
                max[timeAxis] = i;

                //Duplicate datasets to not modify originals when applying masks later
                Dataset dataset1copy = dataset1.duplicate();
                Dataset dataset2copy = dataset2.duplicate();

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


                for (long k = 0; k < radialProfile.Xvalues.length; k++) {
                    for (int l = 0; l < 3; l++) {
                        correlationAccessor.setPosition(new long[]{i,k,l});
                        correlationAccessor.get().setReal(radialProfile.Yvalues[2-l][(int)k]);
                    }
                }

                if(radialProfile.confidence > highestConfidence){
                    highestConfidence = radialProfile.confidence;
                    highestConMean = radialProfile.gaussFit[1];
                    highestConSD = radialProfile.gaussFit[2];
                    highestConFrame = i;
                    highestCCvalue = radialProfile.Yvalues[2][(int)Math.round(radialProfile.gaussFit[1])];
                }

                LinkedHashMap<String, Double> gaussianMap = new LinkedHashMap<String, Double>();
                gaussianMap.put("Mean",  getSigDigits(radialProfile.gaussFit[1]));
                gaussianMap.put("SD",  getSigDigits(radialProfile.gaussFit[2]));
                gaussianMap.put("Gaussian height", getSigDigits(radialProfile.Yvalues[2][(int)Math.round(radialProfile.gaussFit[1])]));
                gaussianMap.put("Confidence",  getSigDigits(radialProfile.confidence));
                listOfGaussianMaps.add(gaussianMap);

                rowNames.add((calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "" + calibratedTime.get().calibratedValue(i) : "Frame " + i));
            }

            if(showIntermediates){
                for (int i = 0; i < 4; i++) {
                    uiService.show(intermediates[i]);
                }
            }

            uiService.show("Heat map of correlation over time between " + dataset1.getName() + " and " + dataset2.getName(), timeCorrelationHeatMap);

            uiService.show("Gaussian fits over time", Tables.wrap(listOfGaussianMaps, rowNames));

            String output = "Highest confidence fit of a gaussian curve to the cross-correlation of: \n\""+
                    dataset1.getName() +
                    "\"\n with \n\"" +
                    dataset2.getName() +
                    "\"\n using the mask \n\"" + (maskAbsent? "No mask selected" : maskDataset.getName()) +
                    "\"\nwas found at " + (calibratedTime.isPresent() && calibratedTime.get().calibratedValue(1) != 0 ? "time " + calibratedTime.get().calibratedValue(highestConFrame) + ", ": "") + "frame " + highestConFrame +
                    ":\n\nMean: " + getSigDigits(highestConMean) +
                    "\nStandard deviation (sigma): " + getSigDigits(highestConSD) +
                    "\nGaussian height: " + getSigDigits(highestCCvalue) +
                    "\nConfidence: " + getSigDigits(highestConfidence) +
                    "\n\n\nThe 3-channel heat map shows the (by channel): \n 1. Gaussian curve for each frame.\n 2. Subtracted correlation for each frame.\n 3. Original correlation for each frame.\n\nFor more details, please see the website: \nhttps://imagej.github.io/Colocalization_by_Cross_Correlation";
            uiService.show("Gauss Fit", output);

            if(highestConfidence < 15){
                uiService.show("Low confidence", "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.");
            }

            if(saveFolder != null) {
                if (!saveFolder.exists() || !saveFolder.canWrite()) {
                    logService.error("Output directory does not exist or does not have write permissions");
                    return;
                }
                try {
                    config.writerSetFailIfOverwriting(false);
                    datasetIOService.save(ContributionOf1, saveFolder.getAbsolutePath() + "\\" + ContributionOf1.getName(), config);
                    datasetIOService.save(ContributionOf2, saveFolder.getAbsolutePath() + "\\" + ContributionOf2.getName(), config);

                    datasetIOService.save(timeCorrelationHeatMap, saveFolder.getAbsolutePath() + "\\" + timeCorrelationHeatMap.getName(), config);

                    ioService.save(Tables.wrap(listOfGaussianMaps, rowNames), saveFolder.getAbsolutePath() + "\\" + "Gaussian fits over time.csv");

                    FileUtils.writeStringToFile(new File(saveFolder.getAbsolutePath() + "\\" + "Results Summary.txt"), output, (Charset) null);

                    if (showIntermediates) {
                        for (Dataset intermediate : intermediates) {
                            datasetIOService.save(intermediate, saveFolder.getAbsolutePath() + "\\" + intermediate.getName(), config);
                        }
                    }

                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        //endregion
    }

    private void colocalizationAnalysis(Img <? extends RealType> img1, Img <? extends RealType> img2, Img <? extends RealType> imgMask, RadialProfiler radialProfiler, final RandomAccessibleInterval <? extends RealType> contribution1, final RandomAccessibleInterval <? extends RealType> contribution2, RandomAccessibleInterval <? extends RealType> [] localIntermediates){
        statusService.showStatus(statusBase + "Applying masks");

        LoopBuilder.setImages(img1, imgMask).multiThreaded().forEachPixel((a,b) -> {if((b.getRealDouble() == 0.0)) {a.setReal(b.getRealDouble());}});
        LoopBuilder.setImages(img2, imgMask).multiThreaded().forEachPixel((a,b) -> {if((b.getRealDouble() == 0.0)) {a.setReal(b.getRealDouble());}});

        statusService.showStatus(statusBase + "Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> oCorr = imgFactory.create(img1);
        Img<FloatType> rCorr = imgFactory.create(img1);
        ExecutorService service = Executors.newCachedThreadPool();

        FFTConvolution conj = new FFTConvolution(img1,img2,service);
        conj.setComputeComplexConjugate(true);
        conj.setOutput(oCorr);
        conj.convolve();

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[0], oCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }


        /** Plot the original correlation in ImageJ as a function of distance
         */

        /**Start creating average correlation of Costes Randomization data. Have to begin this outside the loop to seed
         * avgRandCorr with non-zero data. Data outside the mask is unaltered during the randomization process,
         * which should effectively remove its contribution from the analysis after subtraction. After we initialize
         * it, we can continue to create a average correlation with random data.
         *
         * While working with test data, I noticed that the number of randomizations is not crucial and often a single
         * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
         * data may require more randomization cycles.
         */


        statusService.showStatus(statusBase + "Initializing randomizer");
        CostesRandomizer imageRandomizer = new CostesRandomizer(img1, imgMask);

        Img <RealType> randomizedImage = imageRandomizer.getRandomizedImage();

        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[1], randomizedImage).multiThreaded().forEachPixel((a,b) -> a.setReal(b.getRealFloat()));
        }

        statusService.showStatus(statusBase + "Cycle 1/" + cycles + " - Calculating randomized correlation");
        conj.setImg(randomizedImage);
        conj.setOutput(rCorr);
        conj.convolve();
        Img<FloatType> avgRandCorr = rCorr.copy();

        for (int i = 1; i < cycles; ++i) {
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Randomizing Image");
            randomizedImage = imageRandomizer.getRandomizedImage();
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Calculating randomized correlation");
            conj.setImg(randomizedImage);
            conj.convolve();
            statusService.showStatus(statusBase + "Cycle " + (i+1) + "/" + cycles + " - Averaging randomized correlation");

            ImgMath.compute(ImgMath.div(ImgMath.add(avgRandCorr, rCorr), 2.0)).into(avgRandCorr);
        }

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
            DialogPrompt.Result result = uiService.showDialog("Failed to fit gaussian curve to data, suggesting no correlation between the images.\nSelect OK to continue and show intermediate correlation images. Select cancel to interrupt plugin and show full error message.", DialogPrompt.MessageType.ERROR_MESSAGE, DialogPrompt.OptionType.OK_CANCEL_OPTION);
            if (result == DialogPrompt.Result.CANCEL_OPTION){
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
        ApplyGaussToCorr(subtracted, scale, radialProfiler.Yvalues[2], gaussModifiedCorr);
        if(showIntermediates) {
            LoopBuilder.setImages(localIntermediates[3], gaussModifiedCorr).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
        }
        statusService.showStatus(statusBase + "Determining channel contributions");

        //To get contribution of img1, convolve img2 with the gauss-modified correlation, then multiply with img1
        conj.setComputeComplexConjugate(false);
        conj.setImg(img2);
        conj.setKernel(gaussModifiedCorr);
        conj.setOutput(rCorr);
        conj.convolve();

        LoopBuilder.setImages(contribution1, ImgMath.compute(ImgMath.mul(rCorr, img1)).into(rCorr.copy())).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));

        //To get contribution of img2, correlate img1 with the gauss-modified correlation, then multiply with img2
        conj.setComputeComplexConjugate(true);
        conj.setImg(img1);
        conj.convolve();

        LoopBuilder.setImages(contribution2, ImgMath.compute(ImgMath.mul(rCorr, img2)).into(rCorr.copy())).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));

        service.shutdown();
    }

    private double getSigDigits(double input){ return ((Math.round(input* sigDigits))/ sigDigits);}


    //troubleshooting method for showing images at key points
/*    private void showScaledImg(Img input, String title){
        uiService.show(title, input);
    }*/


    private <T extends RealType> void ApplyGaussToCorr(RandomAccessibleInterval <T> input, double[] scale, Double[] gaussYvalues, RandomAccessibleInterval <T> output){

        double tmax = 0;
        for (Double gaussYvalue : gaussYvalues) {
            if (gaussYvalue > tmax) {
                tmax = gaussYvalue;
            }
        }

        final double max = tmax;

        double scaledValueSq = 0;

        //get image dimensions and center
        int nDims = input.numDimensions();
        if(nDims != scale.length)
            return;
        long [] dims = new long[nDims];
        input.dimensions(dims);

        //Used to recreate bin size to apply gaussian values to image.
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(scale[i],2);
        }
        double binSize = Math.sqrt(scaledValueSq);

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
                    outLooper.get().setReal(looper.get().getRealFloat()*(gaussYvalues[(int)Math.round(Ldistance/binSize)]/max));
                }
            });
        });
    }
}

