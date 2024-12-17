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

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

public class AveragedMask <R extends RealType, F extends FloatType> {

    final Double [] sum = {0.0};
    final Integer [] count = {0};

    public AveragedMask(RandomAccessibleInterval<F> source, RandomAccessibleInterval<R> inputMask){
        sum[0] = 0.0;
        count[0] = 0;
        LoopBuilder.setImages(inputMask, source).multiThreaded().forEachPixel((m, s) -> {
            if ((m.getRealFloat() != 0.0)) {
                synchronized (sum){
                    sum[0] += s.getRealDouble();
                }
                synchronized (count){
                    count[0] += 1;
                }
            }
        });
    }



    public  Img<FloatType> getMaskMeanSubtractedImage(RandomAccessibleInterval<FloatType> source, RandomAccessibleInterval<R> inputMask, ImgFactory<FloatType> floatTypeImgFactory){

        Img<FloatType> returnedImage = floatTypeImgFactory.create(source);

        double mean = getMeanUnderMask();

        LoopBuilder.setImages(source, inputMask, returnedImage).multiThreaded().forEachPixel((s, m, r) -> {
            if ((m.getRealFloat() != 0.0)) {
                r.setReal(s.getRealFloat()-mean);
            }
        });

        return returnedImage;
    }

    public Double getMeanUnderMask(){
        return sum[0]/count[0];
    }

    public int getMaskVoxelCount(){
        return count[0];
    }

}
