

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

import org.scijava.ItemIO;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.DynamicCommand;
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

    @Parameter(label = "No mask (not recommended)?", description = "When checked, performs Costes randomization over the entire image, regardless of what image is selected below.", callback = "maskCallback")
    private boolean maskAbsent;

    @Parameter(label = "Mask: ", description = "The mask over which pixels of image 1 will be randomized. This is important, more details at: imagej.github.io/Colocalization_by_Cross_Correlation", required = false, persist = false)
    private Dataset maskDataset;

    @Parameter(label = "PSF xy(pixel units): ", description = "The width of the point spread function for this image, used for Costes randomization", min = "1")
    private long PSFxy;

    @Parameter(label = "PSF z(pixel units), enter 1 for 2D image:", description = "The depth of the point spread function for this image, used for Costes randomization", min = "1")
    private long PSFz;

    @Parameter(label = "Cycle count: ", description = "The number of Costes randomization cycles to perform. Recommend at least 3, more for sparse signal.", min = "1")
    private long cycles;

    @Parameter(label = "Significant digits: ")
    private int sigDigits;

    @Parameter(label = "Show intermediate images? ", description = "Shows images of numerous steps throughout the algorithm. More details at: imagej.github.io/Colocalization_by_Cross_Correlation")
    private boolean intermediates;

    private Dataset ContributionOf1;
    private Dataset ContributionOf2;


    public Colocalization_by_Cross_Correlation() {
    }

    public void maskCallback(){
        if (maskAbsent) {
            maskDataset = dataset1;
        }
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
        double[] scale = new double[dataset1.numDimensions()];
        for (int i = 0; i < scale.length; i++) {
            scale[i] = dataset1.averageScale(i);
        }

        if(dataset1.getFrames() == 1) {
            RadialProfiler radialProfile = new RadialProfiler(dataset1, scale);
            colocalizationAnalysis(dataset1, dataset2, maskDataset, radialProfile);

            Plot plot = new Plot("Correlation of images","Distance (scaled)", "Relative correlation");
            plot.add("line", radialProfile.Xvalues, radialProfile.Yvalues[0]);
            plot.setColor("red");
            plot.add("line", radialProfile.Xvalues, radialProfile.Yvalues[1]);
            plot.setColor("blue");
            plot.addPoints(radialProfile.Xvalues, radialProfile.Yvalues[2], 0);
            plot.show();
            plot.setLimits(0, radialProfile.gaussFit[1] + (5*radialProfile.gaussFit[2]), 0, plot.getLimits()[3]);
        }

        /**
         * To make time-lapse compatible: write colocalizationAnalysis plot output to the column of a new image one column per time point;
         * needs to be a 3 channel image (Gaussian, oCorr, sCorr)
         * Will have to use the x scale to store time values (read from file if possible?), and the y-axis to store distance
         * (may have to do first iteration out of the loop to determine output image size)
         */
    }

    private <T extends NumericType< T >> void colocalizationAnalysis(Dataset img1, Dataset img2, Dataset imgMask, RadialProfiler radialProfiler){
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
        radialProfiler.calculateProfiles(oCorr, subtracted);

        /**After getting the radial profile, need to fit a gaussian curve to the data, and draw the points to
         * the plot window.
         */

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
        ApplyGaussToCorr(subtracted, scale, radialProfiler.Yvalues[2], gaussModifiedCorr);
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

        uiService.show("Gauss Fit", "Fit a gaussian curve to the cross-correlation of: \n\""+ img1.getName() + "\"\n with \n\"" + img2.getName() + "\"\n using the mask \n\"" + (maskAbsent? "No mask selected" : imgMask.getName()) + "\":\n\nMean: " + Math.round(radialProfiler.gaussFit[1]*significant)/significant + "\nStandard deviation (sigma): " + Math.round(radialProfiler.gaussFit[2]*significant)/significant + "\nConfidence: " + Math.round(radialProfiler.confidence*significant)/significant);

        if(radialProfiler.confidence < 15){
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

    private <T extends RealType> void ApplyGaussToCorr(RandomAccessibleInterval <T> input, double[] scale, double[] gaussYvalues, RandomAccessibleInterval <T> output){

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

