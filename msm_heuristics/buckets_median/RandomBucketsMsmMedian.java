package buckets_median;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import MSMDistances.MSMDistance;

/**
 * random-buckets-msm-median
 */
public class RandomBucketsMsmMedian {
    
    protected final double[][] timeseries;
    protected final int k; 
    protected final double c;
    protected final int amountRandom;
    protected final int bucketsize;
    
    // length of minimum time series
    protected int minDimension = Integer.MAX_VALUE;

    protected double[] randomMean;

    /**
     * @param timeseries
     * @param c 
     * @param amountRandom the amount of random time series
     * @param bucketsize the size of the bucket 
     */
    public RandomBucketsMsmMedian(double[][] timeseries, double c, int amountRandom, int bucketsize) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
        this.amountRandom = amountRandom;
        this.bucketsize = bucketsize;

        for(double[] ts : this.timeseries) {
            this.minDimension = Math.min(minDimension, ts.length);
        }


    }

    /**
     * Selects the best random mean for the time series. 
     * Creates a list of all random time series. 
     */
    public double[][] calc() {
        // list of the random time series
        ArrayList<ArrayList<Double>> randomTs = new ArrayList<>();
        // fill list with random means
        for(int i = 0; i < this.amountRandom; i++) {
            ArrayList<Double> aList = new ArrayList<>();
            int j = getRandomNumber(0, this.bucketsize);
            while(j < minDimension) {
                aList.add(timeseries[getRandomNumber(0, this.k)][j]);
                j = j + getRandomNumber(1, this.bucketsize+1);
            }
            randomTs.add(aList);
        }

        // list of cost of each random time series
        double[] costRandomTs = new double[randomTs.size()];

        MSMDistance msmDistance = new MSMDistance(this.c);

        // iterates through all random means
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

        // selects the time series with the minimum cost
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
        int bucketsize = 2;

        System.out.println();
        double[][] timeseries = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            {2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, 
            {3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
            {4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, 
            {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, 
            {1, 2, 3, 4, 5, 6, 7, 8, 9}, 
            {2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            {3, 4, 5, 6, 7, 8, 9, 10, 11},
            {4, 5, 6, 7, 8, 9, 10, 11, 12}, 
            {5, 6, 7, 8, 9, 10, 11, 12, 13}, 
        };
        int k = timeseries.length;
        System.out.println("Time Series X:");
        for(int i = 0; i < k; i++) {
            System.out.println("x" + (i+1) + ": " + Arrays.toString(timeseries[i]));
        }
        System.out.println();
        System.out.println("Random Buckets Mean: ");
        long startTime = System.nanoTime();
        RandomBucketsMsmMedian rbm = new RandomBucketsMsmMedian(timeseries, c, amountRandom, bucketsize);
        rbm.calc();
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    }

}

