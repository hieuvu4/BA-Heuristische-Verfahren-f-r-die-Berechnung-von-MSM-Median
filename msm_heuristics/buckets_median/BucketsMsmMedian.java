package buckets_median;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import MSMDistances.MSMDistance;
import MSMMean.MsmMean;
import MSMMean.Pair;

/**
 * buckets-msm-median
 */
public class BucketsMsmMedian {
    
    private final double[][] timeseries;
    private final int k;
    private final double c;
    private final int bucketsize;
    private int meanLength = 0;
    private double[] mean;

    private int minDimension = Integer.MAX_VALUE;

    public BucketsMsmMedian(double[][] timeseries, double c, int bucketsize) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
        this.bucketsize = bucketsize;

        if(bucketsize < 1) {
            throw new IllegalArgumentException("Bucketsize must be at least 1");
        }

        for(double[] ts: timeseries) {
            minDimension = Math.min(minDimension, ts.length);
        }       

    }

    // calc msm median for each buckets and adds them together
    public double[][] calc() {
        // List of list of ts with given bucket size
        int i = 0;
        ArrayList<double[]> bucketList = new ArrayList<>();
        while(i < minDimension) {
            // bucket
            double[][] partialTs;
            if(bucketsize > minDimension-i) {
                partialTs = new double[k][minDimension-i];
                for(int j = 0; j < k; j++) {
                    for(int m = 0; m < minDimension-i; m++) {
                        partialTs[j][m] = timeseries[j][i+m];
                    }
                }
                MsmMean msmMean = new MsmMean(partialTs, c);
                Pair<double[], Double> msm = msmMean.calcMean();
                meanLength += msm.getFirst().length;
                bucketList.add(msm.getFirst());
            }
            else {
                partialTs = new double[k][bucketsize];
                for(int j = 0; j < k; j++) {
                    for(int m = 0; m < bucketsize; m++) {
                        partialTs[j][m] = timeseries[j][i+m];
                    }
                }
                MsmMean msmMean = new MsmMean(partialTs, c);
                Pair<double[], Double> msm = msmMean.calcMean();
                meanLength += msm.getFirst().length;
                bucketList.add(msm.getFirst());
            }
            i = i+bucketsize;
        }
        
        i = 0;
        mean = new double[meanLength];
        for(int j = 0; j < bucketList.size(); j++) {
            for(int m = 0; m < bucketList.get(j).length; m++) {
                mean[i] = bucketList.get(j)[m];
                i++;
            }
        }

        double cost = 0.;
        MSMDistance msmDistance = new MSMDistance(c);
        for(int j = 0; j < k; j++) {
            cost += msmDistance.msmDist(timeseries[j], mean);
        }

        double[][] result = new double[2][];
        result[0] = mean;
        double[] costArr = {cost};
        result[1] = costArr;
        return result;
    }


    public static void main(String[] args) {
        double c = 0.1;
        int bucketradius = 4;

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
        BucketsMsmMedian bmm = new BucketsMsmMedian(timeseries, c, bucketradius);
        bmm.calc();
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    
    }

}
