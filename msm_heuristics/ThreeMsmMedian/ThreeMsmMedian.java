package ThreeMsmMedian;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import MSMMean.Pair;
import MSMDistances.MSMDistance;
import MSMMean.MsmMean;

/**
 * 3-msm-median
 */
public class ThreeMsmMedian {

    private final double[][] timeseries;
    private final int k;
    private final double c;

    public ThreeMsmMedian(double[][] timeseries, double c) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
    }


    // calculates median
    public double[][] calc() {
        // check if divisible by 3
        double p = Math.log10(k)/ Math.log10(3);
        if(timeseries.length == 0 || p - (int)p != 0) {
            System.out.println("Set of time series has not a length of n^3");
            return null;
        }
        double[] mean = findMean(timeseries);
        double cost = 0.;
        MSMDistance msmDistance = new MSMDistance(c);
        for(int i = 0; i < k; i++) {
            cost += msmDistance.msmDist(timeseries[i], mean);
        }
        double[][] result = new double[2][];
        result[0] = mean; 
        double[] costArr = {cost};
        result[1] = costArr;
        return result;
    }



    private double[] findMean(double[][] timeseries) {
        if(timeseries.length == 3) {
            MsmMean msmMean = new MsmMean(timeseries, c);
            Pair<double[], Double> msm = msmMean.calcMean();
            return msm.getFirst();
        }
        double[][] ts1 = new double[timeseries.length/3][];
        double[][] ts2 = new double[timeseries.length/3][];
        double[][] ts3 = new double[timeseries.length/3][];
        for(int i = 0; i < timeseries.length/3; i++) {
            ts1[i] = timeseries[i];
            ts2[i] = timeseries[(timeseries.length/3)+i];
            ts3[i] = timeseries[((2 * timeseries.length)/3)+i];
        }

        double[][] lastTs = new double[3][];
        lastTs[0] = findMean(ts1);
        lastTs[1] = findMean(ts2);
        lastTs[2] = findMean(ts3);

        MsmMean msmMean = new MsmMean(lastTs, c);
        Pair<double[], Double> msm = msmMean.calcMean();
        return msm.getFirst();
    }

    public static void main(String[] args) {
        double c = 0.1;

        System.out.println();
        double[][] timeseries = {
            {10, 45, 78, 32, 90},
            {65, 23, 99, 15, 57},
            {81, 3,  41, 74, 68},
            {89, 27, 51, 94, 11},
            {47, 79, 16, 60, 25},
            {63, 14, 92, 7,  83},
            {36, 54, 97, 8,  38},
            {12, 70, 1,  50, 64},
            {42, 22, 85, 58, 19},
            {95, 2,  88, 71, 48},
            {30, 80, 37, 66, 72},
            {4,  31, 20, 86, 26},
            {9,  75, 91, 55, 39},
            {98, 29, 13, 77, 44},
            {67, 5,  84, 17, 53},
            {21, 59, 96, 33, 61},
            {28, 76, 46, 87, 18},
            {40, 82, 34, 56, 69},
            {24, 49, 62, 35, 93},
            {52, 6,  73, 31, 15},
            {3,  67, 42, 90, 11},
            {97, 21, 55, 84, 10},
            {18, 76, 29, 66, 52},
            {40, 1,  83, 48, 26},
            {79, 94, 12, 38, 71},
            {22, 86, 65, 33, 50},
            {7,  57, 30, 99, 73}
        };
        int k = timeseries.length;
        System.out.println("Time Series X:");
        for(int i = 0; i < k; i++) {
            System.out.println("x" + (i+1) + ": " + Arrays.toString(timeseries[i]));
        }
        System.out.println();
        long startTime = System.nanoTime();
        ThreeMsmMedian tMsmMean = new ThreeMsmMedian(timeseries, c);
        double[][] result  = tMsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);

        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    }
    
}
