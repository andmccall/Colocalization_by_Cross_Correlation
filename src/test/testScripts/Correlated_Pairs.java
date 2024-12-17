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
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.DoubleType;


import java.io.IOException;
import java.lang.Math;

public class Correlated_Pairs {

    ImageJ ij;
    protected Dataset paired1,paired2;
    protected String psf1, psf2;


    public Correlated_Pairs(ImageJ inputIJ, String psf1, String psf2){
        ij = inputIJ;
        this.psf1 = psf1;
        this.psf2 = psf2;
    }

    private void getCorrelatedImage(int pairedDistance, int numberOfPoints, int extraPoints, Img<DoubleType> out1, Img<DoubleType> out2){
        long xyDim = out1.dimension(0)-1;
        long zDim = out1.dimension(2)-1;

        RandomAccess<DoubleType> out1accessor = out1.randomAccess();
        RandomAccess<DoubleType> out2accessor = out2.randomAccess();



        for (int point = 0; point < numberOfPoints; point++) {
            long x1 = Math.round(Math.random()*xyDim);
            long y1 = Math.round(Math.random()*xyDim);
            long z1 = Math.round(Math.random()*zDim);

            out1accessor.setPositionAndGet(new long[]{x1, y1, z1}).setReal(out1accessor.get().getRealDouble() + 1.0);

            Double incline = Math.random()*Math.PI;
            Double azimuth = Math.random()*2*Math.PI;

            long x2 = Math.round(pairedDistance*Math.cos(azimuth)*Math.sin(incline))+x1;
            long y2 = Math.round(pairedDistance*Math.sin(azimuth)*Math.sin(incline))+y1;
            long z2 = Math.round(pairedDistance*Math.cos(incline))+z1;

            if(0 <= x2 && x2 < xyDim && 0 <= y2 && y2 < xyDim && 0 <= z2 && z2 <= zDim){
                out2accessor.setPositionAndGet(new long[]{x2, y2, z2}).setReal(out2accessor.get().getRealDouble() + 1.0);
            }
        }

        for (int point = 0; point < extraPoints; point++) {
            long x = Math.round(Math.random()*xyDim);
            long y = Math.round(Math.random()*xyDim);
            long z = Math.round(Math.random()*zDim);

            out1accessor.setPositionAndGet(new long[]{x, y, z}).setReal(out1accessor.get().getRealDouble() + 1.0);

            x = Math.round(Math.random()*xyDim);
            y = Math.round(Math.random()*xyDim);
            z = Math.round(Math.random()*zDim);

            out2accessor.setPositionAndGet(new long[]{x, y, z}).setReal(out2accessor.get().getRealDouble() + 1.0);
        }
    }

    public void generateCorrelatedImages(long xyDim, long zDim, int pairedDistance, int numberOfPoints, int extraPoints) throws IOException {

        Img greenPSF = ij.scifio().datasetIO().open(psf1);
        Img redPSF = ij.scifio().datasetIO().open(psf2);

        Img temp1 = ij.op().create().img(new long[]{xyDim,xyDim,zDim});
        Img temp2 = ij.op().create().img(new long[]{xyDim,xyDim,zDim});

        this.getCorrelatedImage(pairedDistance, numberOfPoints, extraPoints, temp1, temp2);

        paired1 = ij.dataset().create(ij.op().filter().convolve(temp1, greenPSF));
        paired2 = ij.dataset().create(ij.op().filter().convolve(temp2, redPSF));

    }
}
