import java.util.*;

public class MovingAverage {

    private double sum;
    private TreeMap<Double, Double> hashMap;
    private double[] values;

    public MovingAverage(SortedMap<Double,Double> inputMap){
        this.hashMap = new TreeMap<>(inputMap);
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
            output.put(forward.next(), (sum/count));

        }
        return output;
    }

    //Would love to average based on a method like that below, but the computation time is several orders of magnitude greater

/*

        public Map<Double,Double> averagedMap(double range) {
        Map<Double, Double> output = Collections.synchronizedMap(new HashMap<>());
        hashMap.forEach((key,value) -> {
            Double fromKey = key - (range);
            Double toKey = key + (range);
            double newkey = hashMap.subMap(fromKey, toKey).keySet().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            double newvalue = hashMap.subMap(fromKey, toKey).values().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            synchronized (output){output.put(newkey, newvalue);}
        });
        return output;
    }*/

}
