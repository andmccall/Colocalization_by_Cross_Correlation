
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.img.Img;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import java.util.*;
import java.util.Iterator;

public class CostesRandomizer {

    private List<Float> valuesList;

    public CostesRandomizer(Img source, Img inputMask){
        this.setNewImgMask(source, inputMask);
    }

    public void setNewImgMask(Img <? extends RealType> source, Img <? extends RealType> inputMask) {
        valuesList = new ArrayList<>();

        RealType zero = inputMask.getAt(inputMask.minAsLongArray());
        zero.setZero();

        List<Float> syncPositionsList = Collections.synchronizedList(valuesList);

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(inputMask, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor looper = Views.interval(inputMask,chunk).localizingCursor();
                RandomAccess<? extends RealType> sourceAccess = source.randomAccess();
                while(looper.hasNext()){
                    looper.fwd();
                    if(!looper.get().equals(zero.copy())) {
                        Float passer = sourceAccess.setPositionAndGet(looper.positionAsLongArray()).getRealFloat();
                        syncPositionsList.add(passer);
                    }
                }
            });
        });

    }

    public <T extends RealType> Img<T> getRandomizedImage(Img <T> source, Img <T> inputMask){

        Collections.shuffle(valuesList);

        Img randomizedImage = source.copy();

        RealType zero = inputMask.getAt(inputMask.minAsLongArray());
        zero.setZero();

        //List<Float> syncPositionsList = Collections.synchronizedList(valuesList);

        Iterator<Float> iterator = valuesList.iterator();

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(inputMask, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor looper = Views.interval(inputMask,chunk).localizingCursor();
                RandomAccess<? extends RealType> randAccess = randomizedImage.randomAccess();
                while(looper.hasNext()){
                    looper.fwd();
                    if(!looper.get().equals(zero.copy())) {
                        randAccess.setPositionAndGet(looper.positionAsLongArray()).setReal(iterator.next());
                    }
                }
            });
        });

        return randomizedImage;
    }

    public int getMaskVoxelCount(){
        return valuesList.size();
    }

}
