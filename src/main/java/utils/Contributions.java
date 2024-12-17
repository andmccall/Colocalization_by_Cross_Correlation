/*-
 * #%L
 * Scijava plugin for spatial correlation
 * %%
 * Copyright (C) 2019 - 2024 Andrew McCall, University at Buffalo
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
package utils;

import net.imglib2.Cursor;
import net.imglib2.Interval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.algorithm.math.ImgMath;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.view.Views;
import org.apache.commons.math3.analysis.function.Gaussian;

import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Contributions {

    public static void generateGaussianModifiedCCImage(RandomAccessibleInterval<? extends RealType> ccImage, RandomAccessibleInterval <? extends RealType> output, Gaussian gaussian, double [] gaussFitParameters, double [] scale){
        //get image dimensions and center
        int nDims = ccImage.numDimensions();
        if(nDims != scale.length)
            return;
        long [] dims = new long[nDims];
        ccImage.dimensions(dims);

        //obtain center of image
        double[] center = new double[nDims];
        for (int i = 0; i < nDims; i++) {
            center[i] = ((double)dims[i]-1.0)/2;
        }


        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List<Interval> chunks = IntervalChunks.chunkInterval(ccImage, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor<? extends RealType> looper = Views.interval(ccImage,chunk).localizingCursor();
                RandomAccess<? extends RealType> outLooper = output.randomAccess();
                while(looper.hasNext()){
                    looper.fwd();
                    outLooper.setPosition(looper);
                    double LscaledSq = 0;
                    for (int i = 0; i < nDims; ++i) {
                        LscaledSq += Math.pow((looper.getDoublePosition(i)-center[i])*scale[i],2);
                    }
                    double Ldistance = Math.sqrt(LscaledSq);
                    outLooper.get().setReal(looper.get().getRealDouble()*(gaussian.value(Ldistance)/gaussFitParameters[0]));
                    //outLooper.get().setReal(looper.get().getRealDouble()*radialProfile.gaussCurveMap.get(radialProfile.getBD(Ldistance).doubleValue()));
                }
            });
        });
    }

    public static void calculateContributionImages(RandomAccessibleInterval<? extends RealType> img1, RandomAccessibleInterval<? extends RealType> img2, RandomAccessibleInterval <? extends RealType> ccImage, RandomAccessibleInterval<? extends RealType> img1contribution, RandomAccessibleInterval<? extends RealType> img2contribution, ImgFactory<? extends RealType<?>> imgFactory){
        ExecutorService service = Executors.newCachedThreadPool();
        FFTConvolution fdMath = new FFTConvolution(img1, ccImage, img1, imgFactory.imgFactory(new ComplexFloatType()), service);

        calculateContributionImages(img1, img2, ccImage, img1contribution, img2contribution, fdMath);
    }

    public static void calculateContributionImages(RandomAccessibleInterval<? extends RealType> img1, RandomAccessibleInterval<? extends RealType> img2, RandomAccessibleInterval <? extends RealType> ccImage, RandomAccessibleInterval<? extends RealType> img1contribution, RandomAccessibleInterval<? extends RealType> img2contribution, FFTConvolution fdMath){
        fdMath.setComputeComplexConjugate(false);
        fdMath.setImg(img2);
        fdMath.setKernel(ccImage);
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

}

