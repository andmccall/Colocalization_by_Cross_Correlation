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
        oCorrMap = new TreeMap<>();
        sCorrMap = new TreeMap<>();
        gaussCurveMap = new TreeMap<>();

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


        Map<BigDecimal, List<Double>> tempMap2 = new HashMap<>();
        Map<BigDecimal, List<Double>> tempMap = Collections.synchronizedMap(tempMap2);

        //loop through all points, determine distance (scaled)
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
                            tempMap.put(getBD(Ldistance), new ArrayList<>());
                            tempMap.get(getBD(Ldistance)).add(looper.get().getRealDouble());
                        }
                    }
                }
            });
        });

        tempMap.forEach((key,value) -> {
            output.put(key, value.stream().mapToDouble(Double::doubleValue).average().getAsDouble());
        });
        return;
    }


    private double[] CurveFit() {
        WeightedObservedPoints obs = new WeightedObservedPoints();

        sCorrMap.forEach((key, value) -> {
                obs.add(key.doubleValue(), value);
        });

        try {
            return GaussianCurveFitter.create().withMaxIterations(10000).fit(obs.toList());
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

}