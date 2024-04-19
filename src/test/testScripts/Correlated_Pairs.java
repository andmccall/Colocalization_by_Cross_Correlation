import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.real.DoubleType;


import java.io.IOException;
import java.lang.Math;

public class Correlated_Pairs {

    ImageJ ij;

    protected Dataset paired1;
    protected Dataset paired2;

    public Correlated_Pairs(ImageJ inputIJ){
        ij = inputIJ;
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

        Img greenPSF = ij.scifio().datasetIO().open("src/test/resources/green, pre-deconvolved, 10x RL.tif");
        Img redPSF = ij.scifio().datasetIO().open("src/test/resources/red, pre-deconvolved, 10x RL.tif");

        Img temp1 = ij.op().create().img(new long[]{xyDim,xyDim,zDim});
        Img temp2 = ij.op().create().img(new long[]{xyDim,xyDim,zDim});

        this.getCorrelatedImage(pairedDistance, numberOfPoints, extraPoints, temp1, temp2);

        paired1 = ij.dataset().create(ij.op().filter().convolve(temp1, greenPSF));
        paired2 = ij.dataset().create(ij.op().filter().convolve(temp2, redPSF));

    }
}
