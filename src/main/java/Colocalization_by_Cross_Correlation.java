

import ij.*;
import ij.gui.Plot;

import ij.measure.Calibration;
import ij.text.TextWindow;
import net.imagej.Dataset;
import net.imagej.ImgPlus;
import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;

import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.script.math.*;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

import net.imglib2.type.operators.SetOne;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

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

    private String plugTitle = "Colocalization by Cross Correlation";

    @Parameter(label = "Image 1: ")
    protected Dataset dataset1;

    @Parameter(label = "Image 2: ")
    protected Dataset dataset2;

    @Parameter(label = "No mask (not recommended)?", description = "When checked, performs Costes randomization over the entire image")
    protected boolean maskAbsent;

    @Parameter(label = "Mask: ", required = false)
    protected Dataset maskDataset;

    @Parameter(label = "PSF xy(pixel units): ", min = "1")
    protected long PSFxy;

    @Parameter(label = "PSF z(pixel units), enter 1 for 2D image:", min = "1")
    protected long PSFz;

    @Parameter(label = "Cycle count: ", min = "1")
    protected long cycles;

    @Parameter(label = "Show intermediate images? ", description = "Shows ")
    protected boolean intermediates;

    public Colocalization_by_Cross_Correlation() {
    }

    @Override
    public void run(){

        if(dataset1.numDimensions() != maskDataset.numDimensions() && dataset1.numDimensions() != dataset2.numDimensions()){
            IJ.error("Number of image dimensions must be the same");
            return;
        }

        colocalizationAnalysis(dataset1.getImgPlus().getImg(), dataset2.getImgPlus().getImg(), maskDataset.getImgPlus().getImg());
    }

    private <T extends NumericType< T >> void colocalizationAnalysis(Img<? extends NumericType<?>> img1, Img<? extends NumericType<?>> img2, Img<? extends NumericType<?>> imgmask){
        long[] PSF = {PSFxy, PSFxy, PSFz};

        /**wrap input images to img then add those to constructor for FFTConvolutions, then perform an initial
         * correlation of the original data.
         */
//        IJ.showStatus("Wrapping images");

        if(maskAbsent){
            imgmask = img1.copy();
            LoopBuilder.setImages(imgmask).multiThreaded().forEachPixel(SetOne::setOne);
        }


        double[] scale = new double[dataset1.numDimensions()];
        for (int i = 0; i < scale.length; i++) {
            scale[i] = dataset1.averageScale(i);
        }

        IJ.showStatus("Calculating original correlation");

        ImgFactory<FloatType> imgFactory = new ArrayImgFactory<>(new FloatType());
        Img<FloatType> oCorr = imgFactory.create(img1);
        Img<FloatType> rCorr = imgFactory.create(img1);
        ExecutorService service = Executors.newCachedThreadPool();
        FFTConvolution conj = new FFTConvolution(img1,img2,service);
        conj.setComputeComplexConjugate(true);
        conj.setOutput(oCorr);
        conj.convolve();

        if(intermediates) {
            showScaledImg(oCorr, "Original Correlation map", dataset1, scale);
        }

        //Test region to check if correlation of img1 with img2 is mirror of correlation of img2 with img1
        /**
         RandomAccessibleInterval test1 = Views.extendZero(img1);
         RandomAccessible test2 = Views.extendZero(img2);

         conj.setImg(test1);
         conj.setKernel(test2);
         conj.convolve();

         if(intermediates) {
         showScaledImg(oCorr, "Original Correlation map_test1", ip1);
         }

         conj.setImg(test2);
         conj.setKernel(test1);
         conj.convolve();

         if(intermediates) {
         showScaledImg(oCorr, "Original Correlation map_test2", ip1);
         }

         //end of test region
         **/

        /** Plot the original correlation in ImageJ as a function of distance
         */

        Plot oPlot = new Plot("Correlation of images","Distance (scaled)", "Relative correlation");
        IJ.showStatus("Calculating radial profile");
        RadialProfileN(oCorr, scale,oPlot);

        float[][] OdataCurve = oPlot.getDataObjectArrays(0);

        /**Start creating average correlation of Costes Randomization data. Have to begin this outside the loop to seed
         * avgRandCorr with non-zero data. Data outside the mask is unaltered during the randomization process,
         * which should effectively remove its contribution from the analysis after subtraction. After we initialize
         * it, we can continue to create a average correlation with random data.
         *
         * While working with test data, I noticed that the number of randomizations is not crucial and often a single
         * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
         * data may require more randomization cycles.
         */

        IJ.showStatus("Initializing randomizer");
        CostesRandomizer imageRandomizer = new CostesRandomizer(img1, PSF, imgmask);

        Img randomizedImage = imageRandomizer.getRandomizedImage();

        if(intermediates) {
            showScaledImg(randomizedImage, "Initial random image 1", dataset1, scale);
        }
        IJ.showStatus("Cycle 1/" + cycles + " - Calculating randomized correlation");
        conj.setImg(randomizedImage);
        conj.setOutput(rCorr);
        conj.convolve();
        Img<FloatType> avgRandCorr = rCorr.copy();

        for (int i = 1; i < cycles; ++i) {
            IJ.showStatus("Cycle " + (i+1) + "/" + cycles + " - Randomizing Image");
            randomizedImage = imageRandomizer.getRandomizedImage();
            IJ.showStatus("Cycle " + (i+1) + "/" + cycles + " - Calculating randomized correlation");
            conj.setImg(randomizedImage);
            conj.convolve();
            IJ.showStatus("Cycle " + (i+1) + "/" + cycles + " - Averaging randomized correlation data");
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

        IJ.showStatus("Subtracting randomized correlation data");
        Img<FloatType> subtracted;
        try {
            subtracted = new Subtract(oCorr,avgRandCorr).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        if(intermediates) {
            showScaledImg(subtracted, "Subtracted correlation map", dataset1, scale);
        }

        /** Plot the subtracted correlation in the same plot as the original data. Contribution from random elements
         * in the original data should lie close to zero relative to the original data. Real associations will be
         * less affected
         */

        IJ.showStatus("Calculating radial profile");
        oPlot.setColor("red");
        RadialProfileN(subtracted, scale, oPlot);
        float[][] SdataCurve = oPlot.getDataObjectArrays(1);


        oPlot.show();
        oPlot.setLimitsToFit(true);

        /**After getting the radial profile, need to fit a gaussian curve to the data, and draw the points to
         * the plot window.
         */

        IJ.showStatus("Fitting Gaussian");
        double[] Sgaussfit = CurveFit(SdataCurve[0], SdataCurve[1]);

        Gaussian drawCurve = new Gaussian(Sgaussfit[0], Math.abs(Sgaussfit[1]), Sgaussfit[2]);

        float [] gaussYpoints = new float[SdataCurve[0].length];

        for (int i = 0; i < gaussYpoints.length; ++i) {
            gaussYpoints[i] = (float)drawCurve.value(SdataCurve[0][i]);
        }

        oPlot.setColor("blue");
        oPlot.addPoints(SdataCurve[0], gaussYpoints, 0);

        oPlot.setLimits(0, Sgaussfit[1] + (5*Sgaussfit[2]), 0, oPlot.getLimits()[3]);

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

        IJ.showStatus("Calculating confidence");
        double confidence = (areaUnderCurve(SdataCurve[0], SdataCurve[1], Sgaussfit[1], Sgaussfit[2])/areaUnderCurve(OdataCurve[0], OdataCurve[1], Sgaussfit[1], Sgaussfit[2]))*100;

        new TextWindow("Gauss Fit", "Fit a gaussian curve to the correlation of: \n\""+ dataset1.getName() + "\"\n with \n\"" + dataset2.getName() + "\"\n using the mask \n\"" + maskDataset.getName() + "\":\n\nMean: " + Sgaussfit[1] + "\nSigma: " + Sgaussfit[2] + "\nConfidence: " + confidence, 500, 500);


        /** I finally figured it out! To get a representation of the signal from each image that contributed to the
         * cross-correlation after subtraction, I need to do a convolution between the subtracted correlation and
         * img(1?) , then multiply the result of that with the other image.
         *
         * Using rCorr for intermediate steps to avoid generating unnecessary Images
         */

        IJ.showStatus("Determining channel contributions to correlation result");

        //Convolution for img1 and subtracted, then multiply with img2 and overlay
        conj.setComputeComplexConjugate(false);
        conj.setImg(img1);
        conj.setKernel(subtracted);
        conj.setOutput(rCorr);
        conj.convolve();

        Img<FloatType> img2contribution;

        try {
            img2contribution = new Multiply(rCorr, img2).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }
        //This doesn't work, as the order matters from the original correlation, this cannot be accurately reconstructed.
        //Could I flip the correlation data along all the axes to do this?
        showScaledImg(img2contribution, "Contribution of img2", dataset1, scale);

        //Convolution for img2 and subtracted, then multiply with img1 and overlay
        conj.setImg(img2);
        conj.convolve();

        Img<FloatType> img1contribution;

        try {
            img1contribution = new Multiply(rCorr, img1).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        showScaledImg(img1contribution, "Contribution of img1", dataset1, scale);

        //Test for img2 and original, then multiply with img1 and overlay
        conj.setImg(img1);
        conj.setKernel(oCorr);
        conj.convolve();

        Img<FloatType> img3contribution;

        try {
            img3contribution = new Multiply(rCorr, img2).asImage();
        } catch (Exception e) {
            e.printStackTrace();
            return;
        }

        showScaledImg(img3contribution, "Image 2 replicate through reverse correlation", dataset1, scale);

        service.shutdown();

        return;
    }

    private void showScaledImg(RandomAccessibleInterval input, String title, Dataset orig, double[] scaleInfo){
        ImagePlus toShow = ImageJFunctions.wrap(input, title);
        Calibration cal = toShow.getCalibration();

        cal.pixelWidth = scaleInfo[0];
        cal.pixelHeight = scaleInfo[1];
        if(scaleInfo.length > 2)
            cal.pixelDepth = scaleInfo[2];
        toShow.setCalibration(cal);

        toShow.setDimensions( (int)orig.getChannels(), (int)orig.getHeight(), (int)orig.getFrames());
        toShow.setDisplayRange(toShow.getAllStatistics().min, toShow.getAllStatistics().max);
        toShow.show();
        return;
    }

    private <T extends RealType> void RadialProfileN(IterableInterval <T> input, double[] scale, Plot target){

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
        Cursor <T> pointer = input.localizingCursor();

        while(pointer.hasNext()){
            pointer.fwd();
            scaledValueSq = 0;
            for (int i = 0; i < nDims; ++i) {
                scaledValueSq += Math.pow((pointer.getDoublePosition(i)-center[i])*scale[i],2);
            }
            distance = Math.sqrt(scaledValueSq);
            bins[0][(int)Math.round(distance/binSize)] += 1;
            bins[1][(int)Math.round(distance/binSize)] += pointer.get().getRealDouble();
        }



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

}

