import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.List;

public class RadialProfiler {

    //to avoid redundancy, the xvalues get their own array
    public Double[] Xvalues;

    // plotValues[c: 0 = original Correlation, 1 = subtracted correlation, 2 = gaussian fit][binPosition]
    public Double [][]  Yvalues;

    public double[] gaussFit;

    public double confidence;

    private double binSize;

    private long [] dimensions;

    private int nDims;

    private double[] scale;

    public RadialProfiler(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        this.initializeToImageDimensions(input, inputScale);
    }

    static public double getBinSize(RandomAccessibleInterval input, double[] inputScale){
        double scaledValueSq = 0;

        //get image dimensions and center
        int nDims = input.numDimensions();
        double [] scale = inputScale.clone();
        long [] dimensions = new long[nDims];
        input.dimensions(dimensions);

        //Used to create bins of appropriate size for the image.
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(scale[i],2);
        }
        return Math.sqrt(scaledValueSq);
    }

    static public long getNumberOfBins(RandomAccessibleInterval input, double[] inputScale){
        double distance;
        double scaledValueSq = 0;

        //get image dimensions and center
        int nDims = input.numDimensions();
        double [] scale = inputScale.clone();
        long [] dimensions = new long[nDims];
        input.dimensions(dimensions);

        //Used to create bins of appropriate size for the image.
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(scale[i],2);
        }
        double binSize = Math.sqrt(scaledValueSq);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dimensions[i])/2;
        }

        scaledValueSq = 0;
        //calculate the greatest distance from the center, which is also the distance from 0 to the center
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(center[i]*scale[i],2);
        }
        distance = Math.sqrt(scaledValueSq);

        return (int)Math.ceil(distance/binSize)+1;
    }

    //This will need to throw an exception, in case the scale an input dimensions don't match
    public void initializeToImageDimensions(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        double distance;
        double scaledValueSq = 0;

        //get image dimensions and center
        nDims = input.numDimensions();
        if(nDims != inputScale.length)
            throw new Exception("Input image number of dimensions do not match scale number of dimensions");
        scale = inputScale.clone();
        dimensions = new long[nDims];
        input.dimensions(dimensions);

        binSize = getBinSize(input, inputScale);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dimensions[i])/2;
        }

        scaledValueSq = 0;
        //calculate the greatest distance from the center, which is also the distance from 0 to the center
        for (int i = 0; i < nDims; ++i) {
            scaledValueSq += Math.pow(center[i]*scale[i],2);
        }
        distance = Math.sqrt(scaledValueSq);

        Xvalues = new Double[(int)Math.ceil(distance/binSize)+1];

        for (int i = 0; i < Xvalues.length; ++i) {
            Xvalues[i] = binSize * i;
        }
        Yvalues = new Double[3][Xvalues.length];
    }

    public void calculateProfiles(RandomAccessibleInterval origCorrelation, RandomAccessibleInterval subtractedCorrelation){
        calculateSingleProfile(origCorrelation, Yvalues[0]);
        calculateSingleProfile(subtractedCorrelation, Yvalues[1]);
        try{gaussFit = CurveFit();}
        catch (Exception e){
            gaussFit = new double[]{origCorrelation.max(0),origCorrelation.max(0),0};
            throw e;
        }

        Gaussian drawCurve = new Gaussian(gaussFit[0], Math.abs(gaussFit[1]), gaussFit[2]);

        for (int i = 0; i < Xvalues.length; ++i) {
            Yvalues[2][i] = drawCurve.value(Xvalues[i]);
        }

        confidence = (areaUnderCurve(Yvalues[1], gaussFit[1], gaussFit[2])/areaUnderCurve(Yvalues[0], gaussFit[1], gaussFit[2]))*100;
    }

    private <T extends RealType> void calculateSingleProfile(RandomAccessibleInterval <T> input, Double [] output){
        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dimensions[i])/2;
        }

        //bins[0][x] will be count at bin x, bins [1][x] will be integrated density at bin x
        double [][] bins = new double[2][Xvalues.length];

        //loop through all points, determine distance (scaled) and bin

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List<Interval> chunks = IntervalChunks.chunkInterval(input, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor<T> looper = Views.interval(input,chunk).localizingCursor();
                while(looper.hasNext()){
                    looper.fwd();
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i)-center[i])*scale[i],2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    synchronized (bins) {
                        //Have to round half down here for images with odd dimensions that have positive correlation near zero, which can otherwise cause problems with a NaN value for the first bin if half values are rounded up
                        int binPosition = new BigDecimal(Ldistance / binSize).setScale(0, RoundingMode.HALF_DOWN).intValue();
                        bins[0][binPosition] += 1;
                        bins[1][binPosition] += looper.get().getRealDouble();
                        //bins[0][(int) Math.round(Ldistance / binSize)] += 1;
                        //bins[1][(int) Math.round(Ldistance / binSize)] += looper.get().getRealDouble();
                    }
                }
            });
        });

        for (int i = 0; i < Xvalues.length; ++i) {
            output[i] = bins[1][i]/bins[0][i];
        }

        return;
    }


    //The curve is always fit to subtracted, so we use Yvalues[1] throughout
    private double[] CurveFit(){
        int maxLoc = 0;
        double max = 0;

        WeightedObservedPoints obs = new WeightedObservedPoints();

        /**First need to determine the maximum value in order to set the weights for the fitting, and determine its
         * location for instances where the mean is close to zero (in order to mirror the data, this has to be done
         * for a good fit)
         */
        for (int i = 0; i < Xvalues.length; ++i) {
            if (Yvalues[1][i] > max) {
                maxLoc = i;
                max = Yvalues[1][i];
            }
        }

        /** added values are weighted based on the square root of their normalized y-values. The high number of near-zero y-values can
         * mess up the fit
         */

        /**Have to check the possibility of the first & last bin having no values and returning NaN. This is necessary as the
         * gaussian fitter used later will throw an exception with any NaN values.
         */
        for (int i = 0; i < Xvalues.length; ++i) {
            if(!((Double)Yvalues[1][i]).isNaN()) {
                obs.add(Yvalues[1][i] <= 0 ? 0 : Math.sqrt(Yvalues[1][i] / max), Xvalues[i], Yvalues[1][i]);
            }
        }

        /** this next loop adds values below zero that mirror values equidistant from the opposite side of the peak value (max at maxLoc).
         * This is done for fits where the means are near zero, as this data is zero-bounded. Not mirroring the data results
         * in very poor fits for such values. We can't simply mirror across 0 as this will create a double-peak
         * for any data where the peak is near but not at zero.
         * It would be preferable to fit the data using a truncated gaussian fitter, but I could not find any available
         * java class that performs such a fit and my own attempts were unsuccessful.
         */

        for(int i = 1; i < (Xvalues.length - (2*maxLoc)); ++i){
            if(!((Double)Yvalues[1][i + (2 * maxLoc)]).isNaN()) {
                obs.add(Yvalues[1][i + (2 * maxLoc)] <= 0 ? 0 : Math.sqrt(Yvalues[1][i + (2 * maxLoc)] / max), -(Xvalues[i]), Yvalues[1][i + (2 * maxLoc)]);
            }
        }

        try{
            return GaussianCurveFitter.create().withMaxIterations(1000).fit(obs.toList());
        }
        catch(Exception e){
           throw e;
        }

    }

    private double areaUnderCurve(Double[] yvalues, double mean, double sigma){

        double auc = 0;

        for (int i = 0; i < Xvalues.length; ++i) {
            if((mean-(3*sigma)) < Xvalues[i] && Xvalues[i] < (mean+(3*sigma)) ){
                auc += yvalues[i];
            }
        }

        return auc;
    }

}
