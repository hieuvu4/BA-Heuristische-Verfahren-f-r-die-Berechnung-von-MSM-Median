package evolutionary_algorithm;


import java.util.Arrays;
import java.util.TreeSet;
import java.util.concurrent.TimeUnit;

import MSMDistances.MSMDistance;

/**
 * buckets-msm-median-ea
 */
public class BucketMsmMedianEA {

    protected final double[][] timeseries;
    protected final int k;
    protected final double c;
    protected final int generations;
    protected double[] mean;
    protected final int bucketradius;
    
    // all values
    protected final double[] distinctTsValues;

    protected int minDimension = Integer.MAX_VALUE;
    protected int maxDimension = 0;

    public BucketMsmMedianEA(double[][] timeseries, double c, int bucketradius, int generations) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
        this.bucketradius = bucketradius;
        this.generations = generations;

        for(double[] ts: timeseries) {
            this.minDimension = Math.min(this.minDimension, ts.length);
            this.maxDimension = Math.max(this.minDimension, ts.length);
        }
    
        // list of all values
        TreeSet<Double> values = new TreeSet<Double>();
        for (var ts : this.timeseries) {
            for (double v : ts) {
                values.add(v);
            }
        }
        this.distinctTsValues = new double[values.size()];
        int i = 0;
        for (double v : values) {
            this.distinctTsValues[i++] = v;   
        }

        // random mean is generated here
        this.mean = new double[minDimension];
        for(int j = 0; j < minDimension; j++) {
            this.mean[j] = distinctTsValues[getRandomNumber(0, distinctTsValues.length)];
        }
        
        // // random mean is generated here
        // this.mean = new double[minDimension];
        // for(int i = 0; i < this.minDimension; i++) {
        //     this.mean[i] = timeseries[getRandomNumber(0, this.k)][i];
        // }
    }

    /**
     * This is the function where we evolve our mean into a better mean. We iterate 
     * through all generations. In each generation, there is a probability of 1/n for each element 
     * in mean to change into a random number which is in the bucket. The changed mean (newMean) will 
     * be compared with the current mean by calculating the MSM distance with each time series. 
     * If the cost of the changed mean is lower than the cost of the current mean, the current mean 
     * will have the values of the new mean, else nothing happened. 
     */
    public double[][] evolve() {

        System.out.println("Bucket MSM Mean EA: ");
        System.out.println("Time Series Gen. 0: " + Arrays.toString(mean));
        MSMDistance msmDistance = new MSMDistance(c);
        double cost = 0;
        for(double[] ts: timeseries) {
            cost += msmDistance.msmDist(ts, mean);
        }
        double[] newMean = new double[mean.length];

        // go through generations
        for(int i = 0; i < generations; i++) {
            for(int j = 0; j < mean.length; j++) {
                if(getRandomNumber(0, mean.length) == 0) {
                    if(j < bucketradius) {
                        newMean[j] = timeseries[getRandomNumber(0, k)][getRandomNumber(0, bucketradius + 1)];
                    }
                    else if(j >= (mean.length-bucketradius)) {
                        newMean[j] = timeseries[getRandomNumber(0, k)][getRandomNumber(j-bucketradius, mean.length)];
                    }
                    else {
                        newMean[j] = timeseries[getRandomNumber(0, k)][getRandomNumber(j-bucketradius, j+bucketradius)];
                    }
                }
                else {
                    newMean[j] = mean[j];
                }
            }
            double newMeanCost = 0;
            for(double[] ts: timeseries) {
                newMeanCost += msmDistance.msmDist(ts, newMean);
            }
            
            if(newMeanCost < cost) {
                for(int j = 0; j < mean.length; j++) {
                    mean[j] = newMean[j];
                }
                cost = newMeanCost;
            }
            // System.out.println("Time Series Gen. " + i + ": " + Arrays.toString(mean));
            // System.out.println("Cost: " + cost);
        }

        double[][] result = new double[2][];
        result[0] = mean;
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
    protected int getRandomNumber(int min, int max) {
        return (int) ((Math.random() * (max - min)) + min);
    }  
    
    public static void main(String[] args) {
        double c = 0.1;
        int bucketradius = 1;
        int generations = 1_000_000;

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
        long startTime = System.nanoTime();
        BucketMsmMedianEA cmbvea = new BucketMsmMedianEA(timeseries, c, bucketradius, generations);
        cmbvea.evolve();
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    
    }

}
