import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.exception.TooManyIterationsException;
import org.apache.commons.math3.fitting.GaussianCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;

import java.math.BigDecimal;
import java.util.*;


public class RadialProfiler {

    public SortedMap<Double, Double> oCorrMap;
    public SortedMap<Double, Double> sCorrMap;
    public SortedMap<Double, Double> gaussCurveMap;
    //todo: remove gaussCurveMap and just use the Gaussian to generate the data, to free up memory

    public Gaussian gaussian;

    public double[] gaussFitPamameters;

    public double confidence;

    public double rSquared;

    private long[] dimensions;

    private int nDims;

    private double[] scale;

    private int BDscale;

    public RadialProfiler(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        this.initializeToImageDimensions(input, inputScale);
    }

    //This will need to throw an exception, in case the scale an input dimensions don't match
    public void initializeToImageDimensions(RandomAccessibleInterval input, double[] inputScale) throws Exception {
        //get image dimensions and center
        nDims = input.numDimensions();
        if (nDims != inputScale.length)
            throw new Exception("Input image number of dimensions do not match scale number of dimensions");
        scale = inputScale.clone();
        dimensions = new long[nDims];
        input.dimensions(dimensions);

        BDscale = BigDecimal.valueOf(inputScale[0]).scale() + 3;
    }

    public void calculateProfiles(RandomAccessibleInterval origCorrelation, RandomAccessibleInterval subtractedCorrelation) {
        oCorrMap = new TreeMap<>();
        sCorrMap = new TreeMap<>();

        calculateSingleProfile(origCorrelation, oCorrMap);
        calculateSingleProfile(subtractedCorrelation, sCorrMap);
    }

    public <T extends RealType> void calculateSingleProfile(RandomAccessibleInterval<T> input, Map<Double, Double> output) {
        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = (((double) dimensions[i])-1.0) / 2;
        }

        //Map<BigDecimal, Double[]> tempMap2 = new HashMap<>();

        Map<BigDecimal, Double[]> tempMap = Collections.synchronizedMap(new HashMap<>());
        //loop through all points, determine distance (scaled) and bin

        Parallelization.runMultiThreaded(() -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List<Interval> chunks = IntervalChunks.chunkInterval(input, numTasks);

            taskExecutor.forEach(chunks, chunk -> {
                Cursor<T> looper = Views.interval(input, chunk).localizingCursor();
                while (looper.hasNext()) {
                    looper.fwd();
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i) - center[i]) * scale[i], 2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    BigDecimal bdDistance = getBD(Ldistance);
                    synchronized (tempMap) {
                        if (tempMap.containsKey(bdDistance)) {
                            tempMap.get(bdDistance)[0] += looper.get().getRealDouble();
                            tempMap.get(bdDistance)[1] += 1;
                        } else {
                            tempMap.put(bdDistance, new Double[2]);
                            tempMap.get(bdDistance)[0] = looper.get().getRealDouble();
                            tempMap.get(bdDistance)[1] = 1.0;
                        }
                    }
                }
            });
        });

        tempMap.forEach((key,value) -> {
            output.put(key.doubleValue(), (value[0]/value[1]));
        });
    }

    public void fitGaussianCurve(){
        gaussCurveMap = new TreeMap<>();

        gaussFitPamameters = CurveFit();
        gaussian = new Gaussian(gaussFitPamameters[0], Math.abs(gaussFitPamameters[1]), gaussFitPamameters[2]);
        gaussCurveMap.put(gaussFitPamameters[1], gaussian.value(gaussFitPamameters[1]));

        for (Double d : sCorrMap.keySet()) {
            gaussCurveMap.put(d, gaussian.value(d));
        }

        confidence = (areaUnderCurve(sCorrMap, gaussFitPamameters[1], gaussFitPamameters[2]) / areaUnderCurve(oCorrMap, gaussFitPamameters[1], gaussFitPamameters[2]));

        rSquared = getRsquared();

    }

    private double[] CurveFit() {
        double maxLoc = 0;
        double max = 0;
        WeightedObservedPoints obs = new WeightedObservedPoints();

        /*First need to determine the maximum value in order to set the weights for the fitting, and determine its
         * location for instances where the mean is close to zero (in order to mirror the data, this has to be done
         * for a good fit)
         */
        for (Double d : sCorrMap.keySet()) {
            if (sCorrMap.get(d) > max) {
                maxLoc = d;
                max = sCorrMap.get(d);
            }
        }

        if(maxLoc == sCorrMap.firstKey()){
            maxLoc = 0.0;
        }


        /* this next loop adds values below zero that mirror values equidistant from the opposite side of the peak value (max at maxLoc).
         * This is done for fits where the means are near zero, as this data is zero-bounded. Not mirroring the data results
         * in very poor fits for such values. We can't simply mirror across 0 as this will create a double-peak
         * for any data where the peak is near but not at zero.
         * It would be preferable to fit the data using a truncated gaussian fitter, but I could not find any available
         * java class that performs such a fit and my own attempts were unsuccessful.
         */

        if(sCorrMap.firstKey() != 0.0){
            obs.add(0.0, sCorrMap.get(sCorrMap.firstKey()));
        }

        double finalMaxLoc = maxLoc;
        sCorrMap.forEach((key,value) -> {
            obs.add(key, value);
            if (key > 2 * finalMaxLoc) {
                obs.add(((2 * finalMaxLoc) - key), value);
            }
        });

        double [] output = null;
        try{
            output =  GaussianCurveFitter.create().withMaxIterations(100).fit(obs.toList());
        }
        catch(TooManyIterationsException ignored){}


        /*
         * Have to check if the curve was fit to a single noise spike, something that came up quite a bit during initial testing.
         * If not fit to a noise spike, values are returned with no further processing, if it is, the data is averaged
         * with nearest neighbors and another fit is attempted. This usually only needs a single round of averaging.
         *
         * We can use the pixel scale to test for this, as the SD of the spatial correlation should never be less than
         * the pixel size.
         */

        if(output[0] < 0){
            obs.clear();
            double lowest = sCorrMap.values().stream().sorted().findFirst().get();
            if(sCorrMap.firstKey() != 0.0){
                obs.add(0.0, sCorrMap.get(sCorrMap.firstKey()) - lowest);
            }

            sCorrMap.forEach((key,value) -> {
                obs.add(key, value - lowest);
                if (key > 2 * finalMaxLoc) {
                    obs.add(((2 * finalMaxLoc) - key), value - lowest);
                }
            });

            try{
                output =  GaussianCurveFitter.create().withMaxIterations(100).fit(obs.toList());
            }
            catch(TooManyIterationsException ignored){}
        }


        if (output == null || output[2] <= scale[0] || output[1] < 0) {
            //MovingAverage movingAverage = new MovingAverage(sCorrMap);
            for (int windowSize = 1; (output == null || output[2] <= scale[0] || output[1] < 0) && windowSize <= 5; windowSize++) {
                obs.clear();
                SortedMap<Double,Double> averaged = MovingAverage.averagedMap(sCorrMap, windowSize);
                max = 0;
                maxLoc = 0;
                for (Double d : averaged.keySet()) {
                    if (averaged.get(d) > max) {
                        maxLoc = d;
                        max = averaged.get(d);
                    }
                }
                if(maxLoc == averaged.firstKey()){
                    maxLoc = 0.0;
                }
                double finalMaxLoc1 = maxLoc;
                if(averaged.firstKey() != 0.0){
                    obs.add(0.0, averaged.get(averaged.firstKey()));
                }
                averaged.forEach((key, value) -> {
                    obs.add(key, value);
                    if (key > 2 * finalMaxLoc1) {
                        obs.add(((2 * finalMaxLoc1) - key), value);
                    }
                });
                try {
                    output = GaussianCurveFitter.create().withMaxIterations(50).fit(obs.toList());
                }
                catch(TooManyIterationsException ignored){}
            }
        }

        if(output == null|| output[2] <= scale[0] || output[1] < 0){
            gaussFitPamameters = new double[]{0, sCorrMap.lastKey(), sCorrMap.lastKey()};
            confidence = -1;
            rSquared = -1;
            gaussCurveMap = sCorrMap;

            throw new NullPointerException("Could not fit Gaussian curve to data");
        }
        return output;
    }

    private double areaUnderCurve(Map<Double, Double> map, double mean, double sigma) {

        double auc = 0;

        for (Double d : map.keySet()) {
            if ((mean - (3 * sigma)) < d && d < (mean + (3 * sigma))) {
                auc += map.get(d);
            }
        }

        return auc;
    }

    private double getRsquared(){
        final double[] residualsSum = {0};
        final double[] totalSum = {0};
        final double rangeMean = sCorrMap.subMap(gaussFitPamameters[1] - (3 * gaussFitPamameters[2]), gaussFitPamameters[1] + (3 * gaussFitPamameters[2])).values().stream().mapToDouble(Double::doubleValue).average().getAsDouble();

        sCorrMap.subMap(gaussFitPamameters[1] - (3 * gaussFitPamameters[2]),gaussFitPamameters[1] + (3 * gaussFitPamameters[2])).forEach((key,value) -> {
            residualsSum[0] += Math.pow(value - gaussian.value(key), 2);
            totalSum[0] += Math.pow(value - rangeMean, 2);
        });
        return (1-(residualsSum[0]/totalSum[0]));
    }

    public BigDecimal getBD(double input){
        return BigDecimal.valueOf(input).setScale(BDscale, BigDecimal.ROUND_HALF_UP);
    }
}