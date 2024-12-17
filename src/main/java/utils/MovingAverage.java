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

import java.util.*;

public class MovingAverage {

    //Deprecated method, was a little faster, but didn't account for non-linear keys
    /*
    public static SortedMap<Double,Double> averagedMap(SortedMap<Double,Double> hashMap, int windowSize){
        double sum;
        double[] values = hashMap.values().stream().mapToDouble(Double::doubleValue).toArray();

        SortedMap<Double, Double> output = new TreeMap<>();
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

     */

    public static SortedMap<Double,Double> averagedMap(SortedMap<Double,Double> hashMap, Double range) {
        SortedMap<Double, Double> output = Collections.synchronizedSortedMap(new TreeMap<>());
        hashMap.entrySet().parallelStream().forEach((entry) -> {
            Double fromKey = entry.getKey() - range;
            Double toKey = entry.getKey() + range;
            Double newkey = hashMap.subMap(fromKey, toKey).keySet().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            Double newvalue = hashMap.subMap(fromKey, toKey).values().stream().mapToDouble(Double::doubleValue).average().getAsDouble();
            synchronized (output){output.put(newkey, newvalue);}
        });
        return output;
    }
}
