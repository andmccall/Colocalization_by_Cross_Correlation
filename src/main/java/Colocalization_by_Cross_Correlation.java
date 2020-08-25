

import ij.*;
import ij.gui.GenericDialog;
import ij.gui.Plot;

import ij.plugin.PlugIn;
import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.math.Add;
import net.imglib2.algorithm.math.execution.Division;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.script.math.*;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import org.scijava.command.Command;
import org.scijava.plugin.Plugin;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/** An ImageJ co-localization plugin that attempts to find non-random spatial correlations between two images and provide
 * an estimate of their distance and standard deviation. Conceptually similar to Van Steensel's CCF
 *
 * @author Andrew McCall
 */

public class Colocalization_by_Cross_Correlation implements PlugIn{
    //imglib2-script.jar is not included in a FIJI install, needs to be added into Jars folder


    private String plugTitle = "Colocalization by Cross Correlation";

    protected ImagePlus ip1, ip2;
    protected ImagePlus ipmask;
    protected boolean intermediates = false;
    protected long[] PSF = {2,2,2};

    protected long cycles = 3;

    protected String PREF_KEY = "ColocByCorr.";

    public Colocalization_by_Cross_Correlation() {
    }

    public void run(String s){

        if(Dialog()) {
            /**wrap input images to img then add those to constructor for FFTConvolutions, then perform an initial
             * correlation of the original data.
             */
            IJ.showStatus("Wrapping images");
            Img img1 = ImagePlusAdapter.wrap(ip1);
            Img img2 = ImagePlusAdapter.wrap(ip2);
            Img imgmask = ImagePlusAdapter.wrap(ipmask);

            double[] scale = new double[ip1.getNDimensions()];
            scale[0] = ip1.getCalibration().pixelWidth;
            scale[1] = ip1.getCalibration().pixelHeight;
            if (scale.length == 3)
                scale[2] = ip1.getCalibration().pixelDepth;

            for (int i = 0; i < scale.length; ++i) {
                if(scale[i] == 0)
                    scale[i] = 1;
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
                showScaledImg(oCorr, "Original Correlation map", ip1);
            }

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
            IJ.showStatus("Cycle 1/" + cycles + " - Randomizing Image");
            Img randImg = CostesRand(img1, PSF, imgmask);
            if(intermediates) {
                showScaledImg(randImg, "Initial random image 1", ip1);
            }
            IJ.showStatus("Cycle 1/" + cycles + " - Calculating randomized correlation");
            conj.setImg(randImg);
            conj.setOutput(rCorr);
            conj.convolve();
            Img<FloatType> avgRandCorr = rCorr.copy();

            for (int i = 0; i < cycles-1; ++i) {
                IJ.showStatus("Cycle " + (i+2) + "/" + cycles + " - Randomizing Image");
                randImg = CostesRand(img1, PSF, imgmask);
                IJ.showStatus("Cycle " + (i+2) + "/" + cycles + " - Calculating randomized correlation");
                conj.setImg(randImg);
                conj.convolve();
                IJ.showStatus("Cycle " + (i+2) + "/" + cycles + " - Averaging randomized correlation data");
                try {
                    avgRandCorr = new Average(avgRandCorr, rCorr).asImage();
                } catch (Exception e) {
                    e.printStackTrace();
                    return;
                }

            }

            service.shutdown();

            if(intermediates) {
                showScaledImg(avgRandCorr, "Averaged correlation map of Randomized data", ip1);
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
                showScaledImg(subtracted, "Subtracted correlation map", ip1);
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

            new ij.text.TextWindow("Gauss Fit", "Fit a gaussian curve to the correlation of: \n\""+ ip1.getTitle() + "\"\n with \n\"" + ip2.getTitle() + "\"\n using the mask \n\"" + ipmask.getTitle() + "\":\n\nMean: " + Sgaussfit[1] + "\nSigma: " + Sgaussfit[2] + "\nConfidence: " + confidence, 500, 500);


            return;
        }
        return;
    }

    private void showScaledImg(RandomAccessibleInterval input, String title, ImagePlus scaleInfo){
        ImagePlus toShow = ImageJFunctions.wrap(input, title);
        toShow.copyScale(scaleInfo);
        toShow.setDimensions(scaleInfo.getNChannels(), scaleInfo.getNSlices(), scaleInfo.getNFrames());
        toShow.setDisplayRange(toShow.getAllStatistics().min, toShow.getAllStatistics().max);
        toShow.show();

        return;
    }

    private boolean Dialog(){
        String [] winList = WindowManager.getImageTitles();
        String [] maskChoice = new String[winList.length+1];

        maskChoice[0] = "None";
        for (int i=0; i < winList.length; ++i){
            maskChoice[i+1] = winList[i];
        }

        if(winList.length < 2) {
            IJ.error("At least two images and a mask are required.");
            return false;
        }

        PSF[0] = (long)Prefs.get(PREF_KEY + "PSFxy", PSF[0]);
        PSF[2] = (long)Prefs.get(PREF_KEY + "PSFz", PSF[2]);
        cycles = (long)Prefs.get(PREF_KEY + "cycles" ,cycles);
        intermediates = Prefs.get(PREF_KEY + "intermediates", intermediates);


        GenericDialog gd = new GenericDialog(plugTitle);
        gd.addChoice("Image1 (to be randomized):", winList, winList[0]);
        gd.addChoice("Image2:", winList, winList[1]);
        gd.addChoice("Randomization Mask:", maskChoice, maskChoice[0]);
        gd.addNumericField("PSF xy(pixel units):", PSF[0], 0);
        gd.addNumericField("PSF z(pixel units), enter 1 for 2D image:", PSF[2], 0);
        gd.addNumericField("Costes Randomization cycles:", cycles, 0);
        gd.addCheckbox("Show intermediate data", intermediates);

        gd.showDialog();
        if(gd.wasCanceled()) {
            return false;
        }

        ip1 = WindowManager.getImage(gd.getNextChoice());
        ip2 = WindowManager.getImage(gd.getNextChoice());
        String maskString = gd.getNextChoice();
        if(maskString.compareTo("None")==0){
            ipmask = ip1.duplicate();
            for (int i = 0; i < ipmask.getNSlices(); ++i) {
                ipmask.setSlice(i);
                ipmask.getProcessor().min(1);
            }
        }
        else{
            ipmask = WindowManager.getImage(maskString);
        }

        if(ip1.getNDimensions() != ipmask.getNDimensions() && ip1.getNDimensions() != ip2.getNDimensions()){
            IJ.error("Number of image dimensions must be the same");
            return false;
        }


        PSF = new long[3];
        PSF[0] = (long)gd.getNextNumber();
        PSF[1] = PSF[0];
        PSF[2] = (long)gd.getNextNumber();

        cycles = (long)gd.getNextNumber();

        intermediates = gd.getNextBoolean();

        Prefs.set(PREF_KEY + "PSFxy", PSF[0]);
        Prefs.set(PREF_KEY + "PSFz", PSF[2]);
        Prefs.set(PREF_KEY + "cycles" ,cycles);
        Prefs.set(PREF_KEY + "intermediates", intermediates);

        return true;
    }


    private <T extends NumericType<T>> boolean hasAnyNonZero(IterableInterval<T> image){

        /**Returns true if any of the voxels are not zero, returns false otherwise*/
        NumericType<T> zero = image.firstElement().copy();
        zero.setZero();

        Cursor<T> iterator = image.cursor();

        while(iterator.hasNext()){
            /**I have no idea why I have to pass a copy of zero here, but it doesn't work without .copy()*/
            if(!iterator.next().valueEquals(zero.copy())){
                return true;
            }
        }
        return false;
    }

    private <T extends Type <T>> void copyView(IterableInterval source, IterableInterval target){
        /**Copies the source to the target. As the target view came from a copy of the source view, we do not need to
         * worry about image structure (Array vs Cell).
          */
        Cursor<T> sourceC = source.cursor();
        Cursor<T> targetC = target.cursor();
        while(sourceC.hasNext()){
            sourceC.fwd();
            targetC.fwd();
            targetC.get().set(sourceC.get());
        }

    }

    private Img CostesRand(Img orig, long[] size, Img mask){

        //ToDo: Would be better to swap Img to RandomAccessibleInterval, however this class does not have a copy function

        /** This function assumes that orig and mask have the same number of dimensions, and will fail if they do not.
         * As this function extends the data, the image and mask do not need to have the same exact dimensions, though
         * the mask is zero-extended, so everything beyond the original bounds is effectively ignored.

        This function uses Views to move the data around, Views do not contain any of the original data, but effectively
         pointers to the original data, thus by overwriting a block chain of Views that need to be changed, with the data
         of a randomized block-chain, we can effectively and easily randomize the original image
        */

        int ndim = orig.numDimensions();

        Img rImg = orig.copy();

        RandomAccessibleInterval data;
        RandomAccessibleInterval mdata;
        RandomAccessibleInterval target;

        if (ndim == 2){
            /**Convert any 2D image to be a 3D image with a Z-dimension of 1. This makes the remainder of this function
             * much simpler
             */
            data = Views.addDimension(orig, 1, 1);
            mdata = Views.addDimension(mask, 1, 1);
            target = Views.addDimension(rImg, 1, 1);
        }
        else {
            /**There is no method to make a View without modification that I could find, but for simplicity the 3D images
             * need to be in the same state as any 2D image would be at this point, hence the offset by 0
             */
            long[] nullOffset = new long[ndim];
            data = Views.offset(orig, nullOffset);
            mdata = Views.offset(mask, nullOffset);
            target = Views.offset(rImg, nullOffset);
        }


        long[] idim = new long[data.numDimensions()];
        data.dimensions(idim);

        /**Extends the data beyond the original bounds in case a block extends beyond the image, and to allow for the
         * randomization of the starting point of the block cuts.
         * Have to use extendmirrordouble due to possibility of z being one dimension
         */

        RandomAccessible dataView = Views.extendMirrorDouble(data);
        RandomAccessible mdataView = Views.extendZero(mdata);
        RandomAccessible targetView = Views.extendZero(target);

        /**This generates a random negative offset for each dimension, this is to randomize where the blocks cut the data.
         * We offset the view by this random amount, so that the origin is moved slightly.
         */
        long[] rOffset = new long[size.length];
        Random random = new Random();
        for (int i = 0; i < size.length; ++i) {
            rOffset[i] = Math.round(random.nextFloat()*(-size[i]));
        }

        if(ndim ==2)
            rOffset[2] = 0;

        dataView = Views.offset(dataView, rOffset);
        mdataView = Views.offset(mdataView, rOffset);
        targetView = Views.offset(targetView, rOffset);

        //ToDo: attempt to convert the arraylists to sets and parallelize the block add section; may not need to shuffle
        /** Arraylists will store the blocks that are within the mask. origBlocks will be shuffled, then the data copied over to
         * targetBlocks. By using two Arraylists we only have to iterate through the entire image one time.
          */
        ArrayList <FinalInterval> origBlocks = new ArrayList<>();
        ArrayList <FinalInterval> targetBlocks = new ArrayList<>();



        /**Loops through each block of the image and if there are any non-zero voxels within the mask image for that block,
         * then it adds a view of that block to both array lists.
         */

        //Comment out from here to end of for loops to test parallel streams code

        long[] positionMin = new long[mdataView.numDimensions()];
        long[] positionMax = new long[mdataView.numDimensions()];
        RandomAccessibleInterval Looper;

        for(long z = 0; z <= idim[2]+(-rOffset[2])+(idim[2]%size[2]); z += size[2]) {
            positionMin[2] = z;
            positionMax[2] = z + size[2]-1;
            for (long x = 0; x <= idim[0] + (-rOffset[0])+(idim[0]%size[0]); x += size[0]) {
                positionMin[0] = x;
                positionMax[0] = x + size[0]-1;
                for (long y = 0; y <= idim[1] +(-rOffset[1])+(idim[1]%size[1]); y += size[1]) {
                    positionMin[1] = y;
                    positionMax[1] = y + size[1]-1;
                    Looper = Views.interval(mdataView, positionMin, positionMax);
                    if (hasAnyNonZero((IterableInterval)Looper)) {
                        FinalInterval position = new FinalInterval(positionMin,positionMax);
                        origBlocks.add(position);
                        targetBlocks.add(position);
                    }
                }
            }
        }

        /**Shuffle the origBlock array list and copy the data from the shuffled list to the targetBlocks.
         */
        IterableInterval originalSingle;
        IterableInterval targetSingle;
        Collections.shuffle(origBlocks);
        for (int i = 0; i < targetBlocks.size(); ++i) {
            originalSingle = Views.interval(dataView, origBlocks.get(i));
            targetSingle = Views.interval(targetView, targetBlocks.get(i));
            copyView(originalSingle, targetSingle);
        }

        /**At this point, some Costes randomization algorithms will smooth the data, however when I tested this I found that
         * smoothing the data resulted in insufficient subtraction of the correlated images and produced very poor
         * results.
         */

        return rImg;
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

