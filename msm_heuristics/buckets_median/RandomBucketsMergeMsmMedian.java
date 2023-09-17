package buckets_median;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import MSMDistances.MSMDistance;


/**
 * Same as Randombuckets mean. Additionally, there is a variable called mergedRange.
 * If the last element of a bucket is in range of the element of the next bucket, the last element
 * will be overwritten by the new element of the next bucket, else the last element will be not 
 * overwritten and the new element will be added to the list.
 */
public class RandomBucketsMergeMsmMedian extends RandomBucketsMsmMedian{

    private double mergeRange;

    public RandomBucketsMergeMsmMedian(double[][] timeseries, double c, int l, int bucketsize, double mergeRange) {
        super(timeseries, c, l, bucketsize);
        this.mergeRange = mergeRange;
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
                Double d = timeseries[getRandomNumber(0, this.k)][j];
                // if the new value is in merge range, the last element of aList will be overwritten with the new value
                if((aList.size() != 0) && Math.abs(d.doubleValue() - aList.get(aList.size()-1).doubleValue()) <= this.mergeRange) {
                    aList.set(aList.size()-1, d);
                }
                else {
                    aList.add(d);
                }
                j = j + getRandomNumber(1, this.bucketsize);
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

    public static void main(String[] args) {
        double c = 0.1;
        int l = 1_000_000;
        int bucketsize = 2;
        double mergeRange = 2*c;

        System.out.println();
        double[][] timeseries = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 
            {2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, 
            {3, 4, 5, 6, 7, 8, 9, 10, 11, 12},
            {4, 5, 6, 7, 8, 9, 10, 11, 12, 13}, 
            {5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, 
        };
        int k = timeseries.length;
        System.out.println("Time Series X:");
        for(int i = 0; i < k; i++) {
            System.out.println("x" + (i+1) + ": " + Arrays.toString(timeseries[i]));
        }

        System.out.println();
        System.out.println("Random Buckets Mean: ");
        long startTime = System.nanoTime();
        RandomBucketsMergeMsmMedian rbm = new RandomBucketsMergeMsmMedian(timeseries, c, l, bucketsize, mergeRange);
        double[][] result = rbm.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    }
    
}
