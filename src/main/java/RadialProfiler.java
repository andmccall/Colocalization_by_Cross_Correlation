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
import java.util.*;


public class RadialProfiler {


    //to avoid redundancy, the xvalues get their own array
    //public Double[] Xvalues;

    // plotValues[c: 0 = original Correlation, 1 = subtracted correlation, 2 = gaussian fit][binPosition]
    //public Double [][]  Yvalues;

    public SortedMap<BigDecimal, Double> oCorrMap;
    public SortedMap<BigDecimal, Double> sCorrMap;
    public SortedMap<BigDecimal, Double> gaussCurveMap;

    public double[] gaussFit;

    public double confidence;

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

        BDscale = BigDecimal.valueOf(inputScale[0]).scale();
    }

    public void calculateProfiles(RandomAccessibleInterval origCorrelation, RandomAccessibleInterval subtractedCorrelation) {
        oCorrMap = new TreeMap<BigDecimal, Double>();
        sCorrMap = new TreeMap<BigDecimal, Double>();
        gaussCurveMap = new TreeMap<BigDecimal,Double>();

        calculateSingleProfile(origCorrelation, oCorrMap);
        calculateSingleProfile(subtractedCorrelation, sCorrMap);
        try {
            gaussFit = CurveFit();
        } catch (Exception e) {
            gaussFit = new double[]{origCorrelation.max(0), origCorrelation.max(0), 0};
            throw e;
        }
        Gaussian drawCurve = new Gaussian(gaussFit[0], Math.abs(gaussFit[1]), gaussFit[2]);
        gaussCurveMap.put(getBD(gaussFit[1]), drawCurve.value(gaussFit[1]));

        for (BigDecimal d : sCorrMap.keySet()) {
            gaussCurveMap.put(d, drawCurve.value(d.doubleValue()));
        }

        confidence = (areaUnderCurve(sCorrMap, gaussFit[1], gaussFit[2]) / areaUnderCurve(oCorrMap, gaussFit[1], gaussFit[2])) * 100;
    }

    private <T extends RealType> void calculateSingleProfile(RandomAccessibleInterval<T> input, Map<BigDecimal, Double> output) {
        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double) dimensions[i]) / 2;
        }

        //bins[0][x] will be count at bin x, bins [1][x] will be integrated density at bin x
        //double [][] bins = new double[2][Xvalues.length];

        Map<BigDecimal, List<Double>> tempMap2 = new HashMap<>();

        Map<BigDecimal, List<Double>> tempMap = Collections.synchronizedMap(tempMap2);
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
                    synchronized (tempMap) {
                        if (tempMap.containsKey(getBD(Ldistance))) {
                            tempMap.get(getBD(Ldistance)).add(looper.get().getRealDouble());
                        } else {
                            tempMap.put(getBD(Ldistance), new ArrayList<Double>());
                            tempMap.get(getBD(Ldistance)).add(looper.get().getRealDouble());
                        }
                    }
                }
            });
        });

        tempMap.forEach((key,value) -> {
            output.put(key, value.stream().mapToDouble(Double::doubleValue).average().getAsDouble());
            //ij.IJ.log(key + " - " + (value.stream().mapToDouble(Double::doubleValue).sum() / value.size()) + "\n");
        });

   /*     while (iterator.hasNext()) {
            Long key = iterator.next();
            output.put(key, tempMap.get(key).stream().mapToDouble(Double::doubleValue).sum() / tempMap.get(key).size());
            ij.IJ.log(key + " - " + (tempMap.get(key).stream().mapToDouble(Double::doubleValue).sum() / tempMap.get(key).size()) + "\n");
        }*/
        return;
    }


    //The curve is always fit to subtracted, so we use Yvalues[1] throughout
    private double[] CurveFit() {
        double maxLoc = 0;
        double max = 0;
        WeightedObservedPoints obs = new WeightedObservedPoints();

        /**First need to determine the maximum value in order to set the weights for the fitting, and determine its
         * location for instances where the mean is close to zero (in order to mirror the data, this has to be done
         * for a good fit)
         */


        for (BigDecimal d : sCorrMap.keySet()) {
            if (sCorrMap.get(d) > max) {
                maxLoc = d.doubleValue();
                max = sCorrMap.get(d);
            }
        }


        /** this next loop adds values below zero that mirror values equidistant from the opposite side of the peak value (max at maxLoc).
         * This is done for fits where the means are near zero, as this data is zero-bounded. Not mirroring the data results
         * in very poor fits for such values. We can't simply mirror across 0 as this will create a double-peak
         * for any data where the peak is near but not at zero.
         * It would be preferable to fit the data using a truncated gaussian fitter, but I could not find any available
         * java class that performs such a fit and my own attempts were unsuccessful.
         */
        double finalMax = max;
        double finalMaxLoc = maxLoc;
        SortedMap<Double,Double> reflected = new TreeMap<>();
        sCorrMap.keySet().stream().sorted().forEachOrdered((d) -> {
            reflected.put(d.doubleValue(), sCorrMap.get(d));
            if (d.doubleValue() > 2 * finalMaxLoc) {
                reflected.put(((2 * finalMaxLoc) - d.doubleValue()), sCorrMap.get(d));
            }
        });

        /** added values are weighted based on their normalized y-values. The high number of near-zero y-values can
         * mess up the fit
         */
        reflected.forEach((key, value) -> {
                obs.add(key, value);
        });

        try {
            return GaussianCurveFitter.create().fit(obs.toList());
        } catch (Exception e) {
            throw e;
        }
    }

    private double areaUnderCurve(Map<BigDecimal, Double> map, double mean, double sigma) {

        double auc = 0;

        for (BigDecimal d : map.keySet()) {
            if ((mean - (3 * sigma)) < d.doubleValue() && d.doubleValue() < (mean + (3 * sigma))) {
                auc += map.get(d);
            }
        }

        return auc;
    }

    public BigDecimal getBD(double input){
        return BigDecimal.valueOf(input).setScale(BDscale, BigDecimal.ROUND_HALF_UP);
    }

    public static Double getDouble(BigDecimal input){
        return input.doubleValue();
    }
}