package evolutionary_algorithm;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

public class BucketMsmNMedianEA extends BucketMsmMedianEA {

    public BucketMsmNMedianEA(double[][] timeseries, double c, int bucketradius, int generations, int meanLength) {
        super(timeseries, c, bucketradius, generations);

        if (meanLength > ((maxDimension - 1) * this.k + 1)) {
            throw new IllegalArgumentException("The given mean length is too big");
        }

        this.mean = new double[meanLength];
        for(int j = 0; j < meanLength; j++) {
            this.mean[j] = distinctTsValues[getRandomNumber(0, distinctTsValues.length)];
        }
        
    }

    public static void main(String[] args) {
        double c = 0.1;
        int bucketradius = 1;
        int generations = 1_000_000;

        int meanLength = 9;

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
        System.out.println("Change Mean Buckets Value EA: ");
        long startTime = System.nanoTime();
        BucketMsmNMedianEA cmbvea = new BucketMsmNMedianEA(timeseries, c, bucketradius, generations, meanLength);
        cmbvea.evolve();
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    }


    

}
