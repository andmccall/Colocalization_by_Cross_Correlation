import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.SortedMap;


public class MovingAverage {

    private double sum;
    private SortedMap<Double, Double> hashMap;
    private double[] values;

    public MovingAverage(SortedMap<Double,Double> hashMap){
        this.hashMap = hashMap;
        this.values = hashMap.values().stream().mapToDouble(Double::doubleValue).toArray();
    }

    public LinkedHashMap<Double,Double> averagedMap(int windowSize){
        LinkedHashMap<Double, Double> output = new LinkedHashMap<>();
        int position = -1;
        Iterator<Double> forward = hashMap.keySet().iterator();
        while(forward.hasNext()){
            int count = 1;
            ++position;
            sum = values[position];
            for (int i = 1; i <= windowSize; i++) {
                if(position-i >= 0){
                    sum += values[position-i];
                    ++count;
                }
                if(position+i < values.length){
                    sum += values[position+i];
                    ++count;
                }
            }
            output.put(forward.next(), sum/count);
        }
        return output;
    }
}
