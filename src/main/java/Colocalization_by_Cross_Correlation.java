

import ij.gui.Plot;

import io.scif.services.DatasetIOService;
import net.imagej.Dataset;

import net.imagej.ImgPlus;
import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.script.math.*;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.type.operators.SetOne;

import net.imglib2.view.Views;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

@Plugin(type = Command.class, menuPath = "Analyze>Colocalization>Colocalization by Correlation")
public class Colocalization_by_Cross_Correlation implements Command{
    //imglib2-script.jar is not included in a FIJI install, needs to be added into Jars folder

    @Parameter
    private LogService logService;

    @Parameter
    private StatusService statusService;

    @Parameter
    private UIService uiService;

    @Parameter(label = "Image 1: ", description = "This is the image which will be randomized during Costes randomization", persist = false)
    private Dataset dataset1;

    @Parameter(label = "Image 2: ", persist = false)
    private Dataset dataset2;

    @Parameter(label = "No mask (not recommended)?", description = "When checked, performs Costes randomization over the entire image, regardless of what image is selected below.")
    private boolean maskAbsent;

    @Parameter(label = "Mask: ", description = "The mask over which pixels of image 1 will be randomized. This is important, more details at: imagej.github.io/Colocalization_by_Cross_Correlation", required = false, persist = false)
    private Dataset maskDataset;

    @Parameter(label = "PSF xy(pixel units): ", description = "The width of the point spread function for this image, used for Costes randomization", min = "1")
    private long PSFxy;

    @Parameter(label = "PSF z(pixel units), enter 1 for 2D image:", description = "The depth of the point spread function for this image, used for Costes randomization", min = "1")
    private long PSFz;

    @Parameter(label = "Cycle count: ", description = "The number of Costes randomization cycles to perform. Recommend at least 3, more for sparse signal.",min = "1")
    private long cycles;

    @Parameter(label = "Significant digits: ")
    private int sigDigits;

    @Parameter(label = "Show intermediate images? ", description = "Shows images of numerous steps throughout the algorithm. More details at: imagej.github.io/Colocalization_by_Cross_Correlation")
    private boolean intermediates;

    public Colocalization_by_Cross_Correlation() {
    }

    @Override
    public void run(){
        if(dataset1.numDimensions() != dataset2.numDimensions()){
            logService.error("Number of image dimensions must be the same");
            return;
        }
        if(!maskAbsent && dataset1.numDimensions() != maskDataset.numDimensions() ){
            logService.error("Mask dimensions must match image dimensions");
            return;
        }
        if(dataset1.getChannels() > 1 || dataset2.getChannels() > 1 || maskDataset.getChannels() > 1){
            logService.error("Multi-channel images are not supported, requires separate channels");
            return;
        }
        if(dataset1.getFrames() > 1 || dataset2.getFrames() > 1 || maskDataset.getFrames() > 1){
            logService.error("Time-lapse data not yet supported");
            return;
        }
        colocalizationAnalysis(dataset1, dataset2, maskDataset);
    }

    private <T extends NumericType< T >> void colocalizationAnalysis(Dataset img1, Dataset img2, Dataset imgMask){
        long[] PSF = {PSFxy, PSFxy, PSFz};
        double significant = Math.pow(10.0,sigDigits);

        if(maskAbsent){
            imgMask = img1.duplicate();
            LoopBuilder.setImages(imgMask).multiThreaded().forEachPixel(SetOne::setOne);
        }

        double[] scale = new double[img1.numDimensions()];
        for (int i = 0; i < scale.length; i++) {
            scale[i] = img1.averageScale(i);
        }


        statusService.showStatus("Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> oCorr = imgFactory.create(img1);
        Img<FloatType> rCorr = imgFactory.create(img1);
        ExecutorService service = Executors.newCachedThreadPool();
        FFTConvolution conj = new FFTConvolution(img1,img2,service);
        conj.setComputeComplexConjugate(true);
        conj.setOutput(oCorr);
        conj.convolve();

        if(intermediates) {
            showScaledImg(oCorr, "Original cross-correlation result", img1);
        }


        /** Plot the original correlation in ImageJ as a function of distance
         */

        Plot correlationPlots = new Plot("Correlation of images","Distance (scaled)", "Relative correlation");
        statusService.showStatus("Calculating radial profile");
        RadialProfileN(oCorr, scale,correlationPlots);

        float[][] OdataCurve = correlationPlots.getDataObjectArrays(0);

        /**Start creating average correlation of Costes Randomization data. Have to begin this outside the loop to seed
         * avgRandCorr with non-zero data. Data outside the mask is unaltered during the randomization process,
         * which should effectively remove its contribution from the analysis after subtraction. After we initialize
         * it, we can continue to create a average correlation with random data.
         *
         * While working with test data, I noticed that the number of randomizations is not crucial and often a single
         * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
         * data may require more randomization cycles.
         */

        statusService.showStatus("Initializing randomizer");
        CostesRandomizer imageRandomizer = new CostesRandomizer(img1, PSF, imgMask);

        Img randomizedImage = imageRandomizer.getRandomizedImage();

        if(intermediates) {
            showScaledImg(randomizedImage, "Costes randomized image", img1);
        }

        statusService.showStatus("Cycle 1/" + cycles + " - Calculating randomized correlation");
        conj.setImg(randomizedImage);
        conj.setOutput(rCorr);
        conj.convolve();
        Img<FloatType> avgRandCorr = rCorr.copy();

        for (int i = 1; i < cycles; ++i) {
            statusService.showStatus("Cycle " + (i+1) + "/" + cycles + " - Randomizing Image");
            randomizedImage = imageRandomizer.getRandomizedImage();
            statusService.showStatus("Cycle " + (i+1) + "/" + cycles + " - Calculating randomized correlation");
            conj.setImg(randomizedImage);
            conj.convolve();
            statusService.showStatus("Cycle " + (i+1) + "/" + cycles + " - Averaging randomized correlation data");
            try {
                avgRandCorr = new Average(avgRandCorr, rCorr).asImage();
            } catch (Exception e) {
                e.printStackTrace();
                return;
            }

        }

        /**Subtract the random correlation from the original to generate a subtracted correlation map. This is
         * what will be used to evaluate any spatial relations between the two channels.
         */

        statusService.showStatus("Subtracting randomized correlation data");
        Img<FloatType> subtracted;
        try {
            subtracted = new Subtract(oCorr,avgRandCorr).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        if(intermediates) {
            showScaledImg(subtracted, "Subtracted cross-correlation result", img1);
        }

        /** Plot the subtracted correlation in the same plot as the original data. Contribution from random elements
         * in the original data should lie close to zero relative to the original data. Real associations will be
         * less affected
         */

        statusService.showStatus("Calculating radial profile");
        correlationPlots.setColor("red");
        RadialProfileN(subtracted, scale, correlationPlots);
        float[][] SdataCurve = correlationPlots.getDataObjectArrays(1);


        correlationPlots.show();
        correlationPlots.setLimitsToFit(true);

        /**After getting the radial profile, need to fit a gaussian curve to the data, and draw the points to
         * the plot window.
         */

        statusService.showStatus("Fitting Gaussian");
        double[] Sgaussfit = CurveFit(SdataCurve[0], SdataCurve[1]);

        Gaussian drawCurve = new Gaussian(Sgaussfit[0], Math.abs(Sgaussfit[1]), Sgaussfit[2]);

        float [] gaussYpoints = new float[SdataCurve[0].length];

        for (int i = 0; i < gaussYpoints.length; ++i) {
            gaussYpoints[i] = (float)drawCurve.value(SdataCurve[0][i]);
        }

        correlationPlots.setColor("blue");
        correlationPlots.addPoints(SdataCurve[0], gaussYpoints, 0);

        correlationPlots.setLimits(0, Sgaussfit[1] + (5*Sgaussfit[2]), 0, correlationPlots.getLimits()[3]);

        /** Once we have the fit, we need to establish a confidence value in it. This is very important as the
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

        statusService.showStatus("Calculating confidence");
        double confidence = (areaUnderCurve(SdataCurve[0], SdataCurve[1], Sgaussfit[1], Sgaussfit[2])/areaUnderCurve(OdataCurve[0], OdataCurve[1], Sgaussfit[1], Sgaussfit[2]))*100;

        /** I finally figured it out! To get a representation of the signal from each image that contributed to the
         * cross-correlation after subtraction, I need to do a convolution between the subtracted correlation and
         * img(1?) , then multiply the result of that with the other image.
         *
         * Using rCorr for intermediate steps to avoid generating unnecessary Images
         */

        /*
        To get contributions, I can't just use subtracted, as it is equivalent to a lower valued version of the original
        correlation map. I have to modify subtracted by the Gaussian fit results
         */

        Img gaussModifiedCorr = subtracted.copy();
        ApplyGaussToCorr(subtracted, scale, gaussYpoints, gaussModifiedCorr);
        if(intermediates) {
            showScaledImg(gaussModifiedCorr, "Gaussian modified cross-correlation result", img1);
        }
        statusService.showStatus("Determining channel contributions to correlation result");

        //To get contribution of img1, convolve img2 with the gauss-modified correlation, then multiply with img1
        conj.setComputeComplexConjugate(false);
        conj.setImg(img2);
        conj.setKernel(gaussModifiedCorr);
        conj.setOutput(rCorr);
        conj.convolve();

        Img<FloatType> img1contribution;

        try {
            img1contribution = new Multiply(rCorr, img1).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }
        showScaledImg(img1contribution, "Contribution of " + img1.getName(), img1);

        //To get contribution of img2, correlate img1 with the gauss-modified correlation, then multiply with img2
        conj.setComputeComplexConjugate(true);
        conj.setImg(img1);
        conj.convolve();

        Img<FloatType> img2contribution;

        try {
            img2contribution = new Multiply(rCorr, img2).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }
        showScaledImg(img2contribution, "Contribution of " + img2.getName(), img1);

        uiService.show("Gauss Fit", "Fit a gaussian curve to the cross-correlation of: \n\""+ img1.getName() + "\"\n with \n\"" + img2.getName() + "\"\n using the mask \n\"" + (maskAbsent? "No mask selected" : imgMask.getName()) + "\":\n\nMean: " + Math.round(Sgaussfit[1]*significant)/significant + "\nStandard deviation (sigma): " + Math.round(Sgaussfit[2]*significant)/significant + "\nConfidence: " + Math.round(confidence*significant)/significant);

        if(confidence < 15){
            uiService.show("Low confidence", "The confidence value for this correlation is low.\nThis can indicate a lack of significant spatial correlation, or simply that additional pre-processing steps are required.\nFor your best chance at a high confidence value, make sure to:\n\n 1. Use an appropriate mask for your data, and \n\n 2. Perform a background subtraction of your images.\nIdeally the background in the image should be close to zero.");
        }

        service.shutdown();

        return;
    }

    private void showScaledImg(Img input, String title, Dataset orig){
        Dataset toShow = orig.duplicateBlank();
        toShow.setImgPlus(ImgPlus.wrap(input));
        uiService.show(title, input);
        return;
    }

    private <T extends RealType> void RadialProfileN(RandomAccessibleInterval <T> input, double[] scale, Plot target){

        double distance;
        double scaledValueSq = 0;

        //get image dimensions and center
        int nDims = input.numDimensions();
        if(nDims != scale.length)
            return;
        long [] dims = new long[nDims];
        input.dimensions(dims);

        //Used to create bins of appropriate size for the image.
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(scale[i],2);
        }
        double binSize = Math.sqrt(scaledValueSq);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dims[i])/2;
        }

        scaledValueSq = 0;
        //calculate the greatest distance from the center, which is also the distance from 0 to the center
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(center[i]*scale[i],2);
        }
        distance = Math.sqrt(scaledValueSq);

        //bins[0][x] will be count at bin x, bins [1][x] will be integrated density at bin x
        double [][] bins = new double[2][(int)Math.ceil(distance/binSize)+1];

       //loop through all points, determine distance (scaled) and bin

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(input, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor <T> looper = Views.interval(input,chunk).localizingCursor();
                while(looper.hasNext()){
                    looper.fwd();
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i)-center[i])*scale[i],2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    bins[0][(int)Math.round(Ldistance/binSize)] += 1;
                    bins[1][(int)Math.round(Ldistance/binSize)] += looper.get().getRealDouble();
                }
            });
        });

        double[] xvalues = null;
        double[] yvalues = null;

        xvalues = new double[bins[0].length];
        yvalues = new double[bins[0].length];

        for (int i = 0; i < xvalues.length; ++i) {
            yvalues[i] = bins[1][i]/bins[0][i];
            xvalues[i] = binSize * i;
        }

        target.add("line", xvalues,yvalues);

        return;
    }

    private double[] CurveFit(float[] xvals, float[] yvals){
        int maxLoc = 0;
        double max = 0;

        WeightedObservedPoints obs = new WeightedObservedPoints();

        /**First need to determine the maximum value in order to set the weights for the fitting, and determine its
         * location for instances where the mean is close to zero (in order to mirror the data, this has to be done
         * for a good fit)
         */
        for (int i = 0; i < xvals.length; ++i) {
            if (yvals[i] > max) {
                maxLoc = i;
                max = yvals[i];
            }
        }

        /** added values are weighted based on the square root of their normalized y-values. The high number of near-zero y-values can
         * mess up the fit
         */

        /**Have to check the possibility of the first & last bin having no values and returning NaN. This is necessary as the
         * gaussian fitter used later will throw an exception with any NaN values.
         */
        for (int i = 0; i < xvals.length; ++i) {
            if(!((Float)yvals[i]).isNaN()) {
                obs.add(yvals[i] <= 0 ? 0 : Math.sqrt(yvals[i] / max), xvals[i], yvals[i]);
            }
        }

        /** this next loop adds values below zero that mirror values equidistant from the opposite side of the peak value (max at maxLoc).
         * This is done for fits where the means are near zero, as this data is zero-bounded. Not mirroring the data results
         * in very poor fits for such values. We can't simply mirror across 0 as this will create a double-peak
         * for any data where the peak is near but not at zero.
         * It would be preferable to fit the data using a truncated gaussian fitter, but I could not find any available
         * java class that performs such a fit and my own attempts were unsuccessful.
          */

        for(int i = 1; i < (xvals.length - (2*maxLoc)); ++i){
            if(!((Float)yvals[i + (2 * maxLoc)]).isNaN()) {
                obs.add(yvals[i + (2 * maxLoc)] <= 0 ? 0 : Math.sqrt(yvals[i + (2 * maxLoc)] / max), -(xvals[i]), yvals[i + (2 * maxLoc)]);
            }
        }

        return GaussianCurveFitter.create().fit(obs.toList());
    }

    private double areaUnderCurve(float[] xvalues, float[] yvalues, double mean, double sigma){

        double auc = 0;

        for (int i = 0; i < xvalues.length; ++i) {
            if((mean-(3*sigma)) < xvalues[i] && xvalues[i] < (mean+(3*sigma)) ){
                auc += yvalues[i];
            }
        }

        return auc;
    }

    private <T extends RealType> void ApplyGaussToCorr(RandomAccessibleInterval <T> input, double[] scale, float[] gaussYvalues, RandomAccessibleInterval <T> output){

        double tmax = 0;
        for (int i = 0; i < gaussYvalues.length; ++i) {
            if (gaussYvalues[i] > tmax) {
                tmax = gaussYvalues[i];
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
                    outLooper.get().setReal(looper.get().getRealFloat()*(Math.sqrt(gaussYvalues[(int)Math.round(Ldistance/binSize)]/max)));
                }
            });
        });

        return;
    }

}

