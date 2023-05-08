import net.imglib2.img.Img;
import net.imglib2.loops.LoopBuilder;
import net.imglib2.type.numeric.RealType;
import java.util.*;
import java.util.Iterator;
import java.util.stream.IntStream;

public class CostesRandomizer {

    private List<Float> valuesList;

    public CostesRandomizer(Img source, Img inputMask){
        this.setNewImgMask(source, inputMask);
    }

    public void setNewImgMask(Img <? extends RealType> source, Img <? extends RealType> inputMask) {
        valuesList = new ArrayList<>();

        List<Float> syncPositionsList = Collections.synchronizedList(valuesList);

        LoopBuilder.setImages(inputMask, source).multiThreaded().forEachPixel((a, b) -> {
            if ((a.getRealFloat() != 0.0)) {
                syncPositionsList.add(b.getRealFloat());
            }
        });
    }

    public <T extends RealType> Img<T> getRandomizedImage(Img <T> source, Img <T> inputMask){
        Collections.shuffle(valuesList);

        Img<T> randomizedImage = source.copy();
        Iterator<Float> iterator = valuesList.iterator();

        LoopBuilder.setImages(inputMask, randomizedImage).multiThreaded().forEachPixel((a, b) -> {
            if ((a.getRealFloat() != 0.0)) {
                b.setReal(iterator.next());
            }
        });
        return randomizedImage;
    }

    public int getMaskVoxelCount(){
        return valuesList.size();
    }

}
