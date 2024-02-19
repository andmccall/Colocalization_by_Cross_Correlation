import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.math.ImgMath;
import net.imglib2.img.Img;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class CCfunctions {

    private FFTConvolution fdMath;
    private AveragedMask averagedMaskImg1;
    private double maskVolume;
    private double [] scale;

    public <T extends RealType> CCfunctions(Img<T> img1, Img <T> img2, double [] inputScale){
        ExecutorService service = Executors.newCachedThreadPool();
        scale = inputScale.clone();
        maskVolume = 1;
        fdMath = new FFTConvolution(extendImage(img1), img1, extendImage(img2), img2, img1.factory().imgFactory(new ComplexFloatType()), service);
    }

    public <T extends RealType> CCfunctions(Img<T> img1, Img <T> img2, Img <T> mask, double [] inputScale){
        ExecutorService service = Executors.newCachedThreadPool();
        averagedMaskImg1 = new AveragedMask(img1, mask);
        scale = inputScale.clone();
        maskVolume = averagedMaskImg1.getMaskVoxelCount()*getVoxelVolume(inputScale);
        fdMath = new FFTConvolution(extendImage(img1), img1, extendImage(img2), img2, img1.factory().imgFactory(new ComplexFloatType()), service);
    }

    public <F extends FloatType> void calculateCC(Img<F> output){

        //OutOfBoundsFactory zeroBounds = new OutOfBoundsConstantValueFactory<>(0.0);

        //ops.filter().correlate(oCorr, img1, img2, img1.dimensionsAsLongArray(), zeroBounds, zeroBounds);

        fdMath.setComputeComplexConjugate(true);
        fdMath.setOutput(output);
        fdMath.convolve();
        LoopBuilder.setImages(output).multiThreaded().forEachPixel((a) -> a.setReal(a.get()/maskVolume));
    }

    public <T extends RealType, F extends FloatType> void calculateCC(Img<T> img1, Img <T> img2, Img <F> output){
        fdMath.setComputeComplexConjugate(true);
        fdMath.setImg(extendImage(img1), img1);
        fdMath.setKernel(extendImage(img2), img2);
        fdMath.setOutput(output);
        fdMath.convolve();
        LoopBuilder.setImages(output).multiThreaded().forEachPixel((a) -> a.setReal(a.get()/maskVolume));
    }

    /**Start creating average correlation of Pixel Randomization data. The zeroed data outside the mask is unaltered
     * during the randomization process, so that it does not contribute to the result.
     *
     * While working with test data, I noticed that the number of randomizations is not crucial and often a single
     * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
     * data may require more randomization cycles.
     */

    public <T extends RealType> void generateSubtractedCCImage(Img<T> img1, Img <T> img2, Img <T> mask, Img <FloatType> originalCorrelationImg, Img <FloatType> output){
        Img<FloatType> lowFreqComp = originalCorrelationImg.factory().create(img1);

        fdMath.setKernel(extendImage(img2), img2);
        fdMath.setComputeComplexConjugate(true);
        fdMath.setOutput(lowFreqComp);
        fdMath.setImg(extendImage(averagedMaskImg1.getAveragedMask(img1, mask)), img1);
        fdMath.convolve();

        LoopBuilder.setImages(lowFreqComp, originalCorrelationImg, output).multiThreaded().forEachPixel((a, oc, out) -> out.setReal(oc.get()-(a.get()/maskVolume)));

    }

    public <T extends RealType> void generateGaussianModifiedCCImage(RandomAccessibleInterval<T> originalCCImage, RandomAccessibleInterval <T> output, RadialProfiler radialProfile){
        //get image dimensions and center
        int nDims = originalCCImage.numDimensions();
        if(nDims != scale.length)
            return;
        long [] dims = new long[nDims];
        originalCCImage.dimensions(dims);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dims[i]-1.0)/2;
        }


        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List<Interval> chunks = IntervalChunks.chunkInterval(originalCCImage, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor<T> looper = Views.interval(originalCCImage,chunk).localizingCursor();
                RandomAccess<T> outLooper = output.randomAccess();
                while(looper.hasNext()){
                    looper.fwd();
                    outLooper.setPosition(looper);
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i)-center[i])*scale[i],2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    outLooper.get().setReal(looper.get().getRealDouble()*radialProfile.gaussian.value(radialProfile.getBD(Ldistance).doubleValue()));
                    //outLooper.get().setReal(looper.get().getRealDouble()*radialProfile.gaussCurveMap.get(radialProfile.getBD(Ldistance).doubleValue()));
                }
            });
        });
    }

    public <T extends RealType, F extends FloatType> void calculateContributionImages(Img<T> img1, Img <T> img2, RandomAccessibleInterval <F> gaussianModifiedCCImage, RandomAccessibleInterval<T> img1contribution, RandomAccessibleInterval<T> img2contribution){
        fdMath.setComputeComplexConjugate(false);
        fdMath.setImg(img2);
        fdMath.setKernel(gaussianModifiedCCImage);
        fdMath.setOutput(img1contribution);
        fdMath.convolve();

        ImgMath.compute(ImgMath.mul(img1contribution, img1)).into(img1contribution);

        fdMath.setComputeComplexConjugate(true);
        fdMath.setImg(img1);
        fdMath.setOutput(img2contribution);
        fdMath.convolve();
        ImgMath.compute(ImgMath.mul(img2contribution, img2)).into(img2contribution);
        //LoopBuilder.setImages(img2contribution, ImgMath.compute(ImgMath.mul(img2contribution, img2)).into(img2contribution.copy())).multiThreaded().forEachPixel((a,b) -> a.setReal(b.get()));
    }

    public <T extends RealType> Img<T> getRandomizedImage(Img<T> img1, Img<T> imgMask){
        return averagedMaskImg1.getAveragedMask(img1, imgMask);
    }

    //made this to quickly and easily test different extension methods for correlation
    private RandomAccessible extendImage(Img in){
        //return Views.extendMirrorSingle(in); //this is the default method, it causes major issues when there is a flat uniform background (even small numbers) over the whole image with no mask
        //return Views.extendValue(in, ops.stats().median(in).getRealDouble()); //this can cause issues similar to extendMirrorSingle, though slightly less often
        return Views.extendZero(in); //this method seems to be the best for cross-correlation. The original cross-correlation can look terrible with flat background or noise (looks like a pyramid), but this is subtracted out. This method also makes the most intuitive sense, as we don't want to correlate beyond the borders of the image.
    }

    private double getVoxelVolume(double [] scale){
        double volume = 1;
        for (double v : scale) {
            volume *= v;
        }
        return volume;
    }

    //troubleshooting method for showing images at key points
    /*    private void showScaledImg(Img input, String title){
        uiService.show(title, input);
    }*/

}
