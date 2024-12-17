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

import net.imglib2.*;
import net.imglib2.algorithm.fft2.FFTConvolution;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.complex.ComplexFloatType;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.view.Views;
import org.apache.commons.math3.analysis.function.Gaussian;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class CCfunctions <R extends RealType<R>, F extends FloatType> {

    private FFTConvolution fdMath;
    protected AveragedMask averagedMaskImg1;
    private double maskVolume;
    private double [] scale;

    public CCfunctions(RandomAccessibleInterval <FloatType> img1, RandomAccessibleInterval <FloatType> img2, double [] inputScale, ImgFactory<R> imgFactory){
        ExecutorService service = Executors.newCachedThreadPool();
        scale = inputScale.clone();
        maskVolume = 1;
        fdMath = new FFTConvolution(extendImage(img1), img1, extendImage(img2), img2, imgFactory.imgFactory(new ComplexFloatType()), service);
    }

    public CCfunctions(RandomAccessibleInterval<FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<R> mask, double [] inputScale, ImgFactory<R> imgFactory){
        ExecutorService service = Executors.newCachedThreadPool();
        averagedMaskImg1 = new AveragedMask(img1, mask);
        scale = inputScale.clone();
        maskVolume = averagedMaskImg1.getMaskVoxelCount()*getVoxelVolume(inputScale);
        fdMath = new FFTConvolution(extendImage(img1), img1, extendImage(img2), img2, imgFactory.imgFactory(new ComplexFloatType()), service);
    }

    public void calculateCC(RandomAccessibleInterval<F> output){

        //OutOfBoundsFactory zeroBounds = new OutOfBoundsConstantValueFactory<>(0.0);

        //ops.filter().correlate(oCorr, img1, img2, img1.dimensionsAsLongArray(), zeroBounds, zeroBounds);

        fdMath.setComputeComplexConjugate(true);
        fdMath.setOutput(output);
        fdMath.convolve();
        LoopBuilder.setImages(output).multiThreaded().forEachPixel((out) -> out.setReal(out.get()/maskVolume));
    }

    public void calculateCC(RandomAccessibleInterval<FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<FloatType> output){
        fdMath.setComputeComplexConjugate(true);
        fdMath.setImg(extendImage(img1), img1);
        fdMath.setKernel(extendImage(img2), img2);
        fdMath.setOutput(output);
        fdMath.convolve();
        LoopBuilder.setImages(output).multiThreaded().forEachPixel((out) -> out.setReal(out.get()/maskVolume));
    }

    /**Start creating average correlation of Pixel Randomization data. The zeroed data outside the mask is unaltered
     * during the randomization process, so that it does not contribute to the result.
     *
     * While working with test data, I noticed that the number of randomizations is not crucial and often a single
     * randomization results in roughly the same correlation map as 50 randomizations averaged together. Sparse
     * data may require more randomization cycles.
     */

    public void generateSubtractedCCImage(RandomAccessibleInterval<FloatType> img1, RandomAccessibleInterval<FloatType> img2, RandomAccessibleInterval<R> mask, Img <FloatType> output, ImgFactory<FloatType> floatTypeImgFactory){
        Img<FloatType> lowFreqComp = averagedMaskImg1.getMaskMeanSubtractedImage(img1, mask, floatTypeImgFactory);

        fdMath.setKernel(extendImage(img2), img2);
        fdMath.setComputeComplexConjugate(true);
        fdMath.setOutput(output);
        fdMath.setImg(extendImage(lowFreqComp), img1);
        fdMath.convolve();

        LoopBuilder.setImages(output).multiThreaded().forEachPixel((out) -> out.setReal(out.get()/maskVolume));
    }

    public void generateGaussianModifiedCCImage(RandomAccessibleInterval<R> ccImage, RandomAccessibleInterval <R> output, Gaussian gaussian, double [] gaussFitParameters){
        Contributions.generateGaussianModifiedCCImage(ccImage, output, gaussian, gaussFitParameters, scale);
    }

    public void calculateContributionImages(RandomAccessibleInterval<? extends RealType> img1, RandomAccessibleInterval<? extends RealType> img2, RandomAccessibleInterval <? extends FloatType> ccImage, RandomAccessibleInterval<? extends RealType> img1contribution, RandomAccessibleInterval<? extends RealType> img2contribution){
        Contributions.calculateContributionImages(img1, img2, ccImage, img1contribution, img2contribution, fdMath);
    }


                                                //made this to quickly and easily test different extension methods for correlation
    private RandomAccessible extendImage(RandomAccessibleInterval<FloatType> in){
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
