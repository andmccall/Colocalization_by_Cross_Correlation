
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.view.Views;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;

public class CostesRandomizer {

    private long[] imageDimensions;
    private Img source;
    private List<long[]> PositionsList, RandomizedList;

    public CostesRandomizer(Img input, Img inputMask){
        this.setNewImg(input, inputMask);
    }

    public <T extends NumericType< T >> void setNewImg(Img <T> input, Img <T> inputMask) {

        imageDimensions = new long[inputMask.numDimensions()];
        inputMask.dimensions(imageDimensions);

        source = input.copy();

        PositionsList = new ArrayList<>();

        NumericType<T> zero = inputMask.getAt(inputMask.minAsLongArray());
        zero.setZero();

        List<long[]> syncPositionsList = Collections.synchronizedList(PositionsList);

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(inputMask, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor looper = Views.interval(inputMask,chunk).localizingCursor();
                while(looper.hasNext()){
                    looper.fwd();
                    if(!looper.get().equals(zero.copy())) {
                        syncPositionsList.add(looper.positionAsLongArray());
                    }
                }
            });
        });

        RandomizedList = new ArrayList<>();
        RandomizedList.addAll(PositionsList);
    }

    public <T extends Type<T>> Img getRandomizedImage(){
        Img randomizedImage = source.copy();

        Collections.shuffle(RandomizedList);

        IntStream.range(0, PositionsList.size()).parallel().forEach(i ->{
            RandomAccess<T> sourcePoint = source.randomAccess(), targetPoint = randomizedImage.randomAccess();
            sourcePoint.setPosition(PositionsList.get(i));
            targetPoint.setPosition(RandomizedList.get(i));
            targetPoint.get().set(sourcePoint.get());
        });

        return randomizedImage;
    }

    public int getMaskVoxelCount(){
        return PositionsList.size();
    }

}
