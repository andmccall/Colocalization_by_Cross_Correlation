
import net.imglib2.*;
import net.imglib2.img.Img;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.loops.IntervalChunks;
import net.imglib2.parallel.Parallelization;
import net.imglib2.parallel.TaskExecutor;
import net.imglib2.realtransform.RealViews;
import net.imglib2.realtransform.Scale2D;
import net.imglib2.realtransform.Scale3D;
import net.imglib2.type.Type;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.view.ExtendedRandomAccessibleInterval;
import net.imglib2.view.Views;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CostesRandomizer {

    protected long[] imageDimensions, PSFsize;
    protected Img source;
    protected List<long[]> PositionsList, RandomizedList;

    public CostesRandomizer(Img input, long[] inputPSFsize, Img inputMask){
        this.setNewImg(input, inputPSFsize, inputMask);
    }

    public <T extends NumericType< T >> void setNewImg(Img <T> input, long[] inputPSFsize, Img <T> inputMask) {

        PSFsize = inputPSFsize.clone();

        imageDimensions = new long[inputMask.numDimensions()];
        inputMask.dimensions(imageDimensions);

        source = input.copy();

        long [] bufferedDimensionsMax = new long[inputMask.numDimensions()];
        long [] bufferedDimensionsMin = new long[inputMask.numDimensions()];

        for (int i = 0; i < inputMask.numDimensions(); i++) {
            bufferedDimensionsMin[i] = 0;
            bufferedDimensionsMax[i] = imageDimensions[i]-1 + (imageDimensions[i]%PSFsize[i] == 0? 0 : PSFsize[i] - (imageDimensions[i]%PSFsize[i]));
        }

        NLinearInterpolatorFactory factory = new NLinearInterpolatorFactory<T>();

        RealRandomAccessible <T> extendedMask = Views.interpolate(Views.interval(Views.expandBorder(inputMask, PSFsize), bufferedDimensionsMin, bufferedDimensionsMax), factory);

        RandomAccessible <T> scalar;

       if (inputMask.numDimensions() == 2){
           scalar = RealViews.affine(extendedMask, new Scale2D(1.0/PSFsize[0], 1.0/PSFsize[1]));
       }
        else if(inputMask.numDimensions() == 3) {
           scalar = RealViews.affine(extendedMask, new Scale3D(1.0 / PSFsize[0], 1.0 / PSFsize[1], 1.0 / PSFsize[2]));
       }
        else{
            return;
       }
        for (int i = 0; i < inputMask.numDimensions(); i++) {
            bufferedDimensionsMax[i] = bufferedDimensionsMax[i]/PSFsize[i];
        }

        RandomAccessibleInterval <T> finalScalar = Views.interval(scalar, bufferedDimensionsMin, bufferedDimensionsMax);

        PositionsList = new ArrayList<>();

        NumericType<T> zero = finalScalar.getAt(finalScalar.minAsLongArray());
        zero.setZero();

        List<long[]> syncPositionsList = Collections.synchronizedList(PositionsList);

        Parallelization.runMultiThreaded( () -> {
            TaskExecutor taskExecutor = Parallelization.getTaskExecutor();
            int numTasks = taskExecutor.suggestNumberOfTasks();
            List< Interval > chunks = IntervalChunks.chunkInterval(finalScalar, numTasks );

            taskExecutor.forEach(chunks, chunk ->{
                Cursor looper = Views.interval(finalScalar,chunk).localizingCursor();
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




    public Img getRandomizedImage(){
        Img randomizedImage = source.copy();

        //Extends the data beyond the original bounds in case a block extends beyond the image
        ExtendedRandomAccessibleInterval data = Views.extendMirrorSingle(source);
        ExtendedRandomAccessibleInterval target = Views.extendMirrorSingle(randomizedImage);

        //Shuffle one of the lists before pairing the two lists to each other
        Collections.shuffle(RandomizedList);

        //Collection map needed for multithreaded processing
        Map<long[], long[]> map = IntStream.range(0, PositionsList.size()).boxed().collect(Collectors.toMap(PositionsList::get, RandomizedList::get));

        //Copies the data from each key position as chunk of the map to each value position as chunk
       map.entrySet().parallelStream().forEach((i) -> {
           long [] sourcemin = new long[i.getKey().length],  sourcemax = new long[i.getKey().length];
           long [] targetmin = new long[i.getValue().length],  targetmax = new long[i.getValue().length];

           //Have to rescale the values in PositionsList as they were generated from a scaled mask

           rescale(i.getKey(), sourcemin, sourcemax);
           rescale(i.getValue(), targetmin, targetmax);

           IterableInterval originalSingle = Views.interval(data, sourcemin, sourcemax);
           IterableInterval targetSingle = Views.interval(target, targetmin, targetmax);
           copyView(originalSingle, targetSingle);

       });


        /**At this point, some Costes randomization algorithms will smooth the data, however when I tested this I found that
         * smoothing the data resulted in insufficient subtraction of the correlated images and produced very poor, inaccurate
         * results.
         */

        return randomizedImage;
    }

    //returns the position in as a rescaled interval defined by min and max, and of size PSFsize
    private void rescale(long[] in, long[] min, long[] max){
        for (int i = 0; i < in.length; i++) {
            min[i] = in[i]*PSFsize[i];
            max[i] = (in[i]*PSFsize[i]) + PSFsize[i]-1;
        }
    }

    private <T extends Type<T>> void copyView(IterableInterval source, IterableInterval target){
        /**Copies the source to the target. As the target view came from a copy of the source view, we do not need to
         * worry about image structure (Array vs Cell).
         */
        Cursor<T> sourceC = source.cursor();
        Cursor<T> targetC = target.cursor();
        while(sourceC.hasNext()){
            sourceC.fwd();
            targetC.fwd();
            targetC.get().set(sourceC.get());
        }

    }
}
