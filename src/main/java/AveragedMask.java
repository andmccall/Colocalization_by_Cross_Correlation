import net.imglib2.img.Img;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import java.util.*;
import java.util.Iterator;
import java.util.stream.IntStream;

public class AveragedMask {

    final Double [] sum = {0.0};
    final Integer [] count = {0};

    public <T extends RealType> AveragedMask(Img <T> source, Img <T> inputMask){
        sum[0] = 0.0;
        count[0] = 0;
        LoopBuilder.setImages(inputMask, source).multiThreaded().forEachPixel((a, b) -> {
            if ((a.getRealFloat() != 0.0)) {
                synchronized (sum){
                    sum[0] += b.getRealDouble();
                }
                synchronized (count){
                    count[0] += 1;
                }
            }
        });
    }

    public <T extends RealType> Img<T> getAveragedMask(Img <T> source, Img <T> inputMask){

        Img<T> returnedImage = source.copy();

        double mean = sum[0]/count[0];

        LoopBuilder.setImages(inputMask, returnedImage).multiThreaded().forEachPixel((a, b) -> {
            if ((a.getRealFloat() != 0.0)) {
                b.setReal(mean);
            }
        });

        return returnedImage;
    }

    public int getMaskVoxelCount(){
        return count[0];
    }

}
