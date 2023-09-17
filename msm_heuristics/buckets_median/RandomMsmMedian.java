package buckets_median;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;

// import MSMDistance class by Jana Holznigenkemper
import MSMDistances.MSMDistance;

public class RandomMsmMedian {
    
    private final double[][] timeseries;
    private final int k;
    private final double c;
    protected final int amountRandom;
    protected double[] distinctTsValues;

    private int minDimension = Integer.MAX_VALUE;

    public RandomMsmMedian(double[][] timeseries, double c, int amountRandom) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
        this.amountRandom = amountRandom;

        for(double[] ts: timeseries) {
            minDimension = Math.min(minDimension, ts.length);
        }       

        TreeSet<Double> values = new TreeSet<Double>();
        for (double[] ts : timeseries) {
            for (double v : ts) {
                values.add(v);
            }
        }
        this.distinctTsValues = new double[values.size()];
        int i = 0;
        for (double v : values)
            this.distinctTsValues[i++] = v;
        // System.out.println(Arrays.toString(distinctTsValues));

    }

    public double[][] calc() {
        // List of list of ts with given bucket size
        ArrayList<ArrayList<Double>> randomTs = new ArrayList<>();
        for(int i = 0; i < this.amountRandom; i++) {
            ArrayList<Double> aList = new ArrayList<>();
            for(int j = 0; j < minDimension; j++) {
                
                aList.add(distinctTsValues[getRandomNumber(0, distinctTsValues.length)]);
            }
            randomTs.add(aList);
        }

        // list of cost of each random time series
        double[] costRandomTs = new double[randomTs.size()];

        MSMDistance msmDistance = new MSMDistance(this.c);

        for(int i = 0; i < this.amountRandom; i++) {
            double sum = 0.;
            double[] currentRandomTs = randomTs.get(i).stream().mapToDouble(Double::doubleValue).toArray();
            // System.out.println(Arrays.toString(currentRandomTs));

            // calculates MSM Distance between all given time series and the current random mean
            for(int j = 0; j < this.k; j++) {            
                sum += msmDistance.msmDist(this.timeseries[j], currentRandomTs);
            }
            costRandomTs[i] = sum;
            // System.out.println(sum);
        }

        double cost = Double.MAX_VALUE;
        int index = 0;
        for(int i = 0; i < costRandomTs.length; i++) {
            if(costRandomTs[i] < cost) {
                cost = costRandomTs[i];
                index = i;
            }
        }

        double[][] result = new double[2][];
        result[0] = randomTs.get(index).stream().mapToDouble(Double::doubleValue).toArray();
        double[] costArr = {cost};
        result[1] = costArr;
        
        return result;
    }



    /**
     * Returns a random integer between @param min (inclusive) and @param max (exclusive).
     * @param min The minimum random integer.
     * @param max The upper bound of the random numbers.
     * @return a random integer.
     */
    protected static int getRandomNumber(int min, int max) {
        return (int) ((Math.random() * (max - min)) + min);
    }   

    public static void main(String[] args) {
        double c = 0.1;
        int amountRandom = 1_000_000;

        System.out.println();
        double[][] timeseries = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            {2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, 
            {3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
            {4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, 
            {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, 
            // {1, 2, 3, 4, 5, 6, 7, 8, 9}, 
            // {2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            // {3, 4, 5, 6, 7, 8, 9, 10, 11},
            // {4, 5, 6, 7, 8, 9, 10, 11, 12}, 
            // {5, 6, 7, 8, 9, 10, 11, 12, 13}, 
        };
        int k = timeseries.length;
        System.out.println("Time Series X:");
        for(int i = 0; i < k; i++) {
            System.out.println("x" + (i+1) + ": " + Arrays.toString(timeseries[i]));
        }
        System.out.println();
        long startTime = System.nanoTime();
        RandomMsmMedian bmm = new RandomMsmMedian(timeseries, c, amountRandom);
        double[][] result = bmm.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println(result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    
    }

}
