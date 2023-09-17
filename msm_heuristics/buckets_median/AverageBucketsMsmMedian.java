package buckets_median;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import MSMDistances.MSMDistance;


public class AverageBucketsMsmMedian {
    
    private final double[][] timeseries;
    private final int k;
    private double c;
    private final int amountRandom;
    private final int bucketSize;
    private final int bucketParts;
    private final int iteration;
    private final double[][][] partialTimeseries;
    private double[] minDistValues;
    private double[][] minDistMean;

    public double[] mean;

    
    // an array which represents the values of each bucket
    private double[][] values;

    private int minDimension = Integer.MAX_VALUE;

    public AverageBucketsMsmMedian(double[][] timeseries, double c, int amountRandom, int bucketSize, int iteration) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;
        this.amountRandom = amountRandom;
        this.bucketSize = bucketSize;
        this.iteration = iteration;
        for(double[] ts : timeseries) {
            minDimension = Math.min(minDimension, ts.length);
        }
        mean = new double[minDimension];

        bucketParts = (minDimension % this.bucketSize == 0) ?  minDimension/this.bucketSize : minDimension/this.bucketSize + 1;

        minDistValues = new double[bucketParts];
        minDistMean = new double[bucketParts][];
        for(int i = 0; i < minDistValues.length; i++) {
            minDistValues[i] = Double.MAX_VALUE;
        }
        values = new double[bucketParts][];

        double[][] rotatedTs = rotateArrayCW(timeseries);
        int currentTsIndex = 0;
        // initialize values
        for(int i = 0; i < values.length; i++) {
                int currentBucketSize = (minDimension - currentTsIndex < this.bucketSize) ? minDimension - currentTsIndex : this.bucketSize;
                double[][] currentBucket = new double[currentBucketSize][k];
            for(int j = 0; j < currentBucketSize; j++) {
                currentBucket[j] = rotatedTs[currentTsIndex];
                currentTsIndex++;
            }

            TreeSet<Double> valuesOfBuckets = new TreeSet<Double>();
            for (var ts : currentBucket) {
                for (double v : ts) {
                    valuesOfBuckets.add(v);
                }
            }

            values[i] = new double[valuesOfBuckets.size()];
            int j = 0;
            for(double v : valuesOfBuckets) {
                values[i][j++] = v;
            }
        }

        // for(double[] val : values) {
        //     System.out.println(Arrays.toString(val));
        // }

        // initialize partial time series 
        currentTsIndex = 0;
        partialTimeseries = new double[bucketParts][][];
        for(int i = 0; i < bucketParts; i++) {
            partialTimeseries[i] = new double[k][];
            int currentBucketSize = (minDimension - currentTsIndex  < this.bucketSize) ? minDimension - currentTsIndex : this.bucketSize;
            for(int j = 0; j < k; j++) {
                partialTimeseries[i][j] = new double[currentBucketSize];
                for(int m = 0; m < currentBucketSize; m++) {
                    partialTimeseries[i][j][m] = timeseries[j][currentTsIndex];
                    currentTsIndex++;
                }
                currentTsIndex -= currentBucketSize;
            }
            currentTsIndex += currentBucketSize;
        }

        // for(int i = 0; i < partialTimeseries.length; i++) {
        //     for(int j = 0; j < partialTimeseries[i].length; j++) {
        //         System.out.println(Arrays.toString(partialTimeseries[i][j]));
        //     }
        //     System.out.println("--");
        // }

    }

    public double[][] calc() {
        for(int i = 0; i < bucketParts; i++) {
            generateRandomBucketMean(partialTimeseries[i], i);
        }
        // for(int i = 0; i < minDistMean.length; i++) {
        //     System.out.println(Arrays.toString(minDistMean[i]));
        //     System.out.println(minDistValues[i]);
        // }

        for(int i = 0; i < iteration; i++) {
            double average = Arrays.stream(minDistValues).average().orElse(Double.NaN);
            for(int j = 0; j < bucketParts; j++) {
                if(minDistValues[j] > average) {
                    generateRandomBucketMean(partialTimeseries[j], j);
                }
            }
        }

        int meanIndex = 0;
        for(int i = 0; i < minDistMean.length; i++) {
            for(int j = 0; j < minDistMean[i].length; j++) {
                mean[meanIndex] = minDistMean[i][j];
                meanIndex++;
            }
        }

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
    

    private void generateRandomBucketMean(double[][] partialTimeseries, int bucketIndex) {
        ArrayList<double[]> randomPartMeans = new ArrayList<>();
        for(int i = 0; i < amountRandom; i++) {
            double[] randomMean = new double[partialTimeseries[0].length];
            for(int j = 0; j < randomMean.length; j++) {
                randomMean[j] = values[bucketIndex][getRandomNumber(0, values[bucketIndex].length)];
            }
            randomPartMeans.add(randomMean);
        }

        double currentDistValue = 0.;
        int minIndex = 0;
        for(int i = 0; i < randomPartMeans.size(); i++) {
            MSMDistance msmDistance = new MSMDistance(c);
            for(int j = 0; j < k; j++) {
                currentDistValue += msmDistance.msmDist(partialTimeseries[j], randomPartMeans.get(i));
            }
            if(currentDistValue < minDistValues[bucketIndex]) {
                minDistValues[bucketIndex] = currentDistValue;
                minIndex = i;
                minDistMean[bucketIndex] = randomPartMeans.get(minIndex);
            }
        }
        
    }
    

    /**
     * This function rotates an 2D double array clockwise
     * @param arr the 2D double array which should be rotated
     * @return a clockwise rotated 2D double array
     */
    private static double[][] rotateArrayCW(double[][] arr) {
        final int rows = arr.length;
        final int cols = arr[0].length;

        final double[][] rotArr = new double[cols][rows];

        for ( int r = 0; r < rows; r++ ) {
            for ( int c = 0; c < cols; c++ ) {
                rotArr[c][rows-1-r] = arr[r][c];
            }
        }
        return rotArr;
    }

    /**
     * Returns a random integer between @param min (inclusive) and @param max (exclusive)
     * @param min The minimum random integer
     * @param max The upper bound of the random numbers
     * @return a random integer
     */
    private static int getRandomNumber(int min, int max) {
        return (int) ((Math.random() * (max - min)) + min);
    }  

    public static void main(String[] args) {

        double c = 0.1;
        int amountRandom = 1_000_000;
        int bucketsize = 2;
        int iteration = 10;

        double[][] timeseries = {
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, 
            {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, 
            {3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13},
            {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14}, 
            {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, 
        };
        
        AverageBucketsMsmMedian abm = new AverageBucketsMsmMedian(timeseries, c, amountRandom, bucketsize, iteration);
        abm.calc();
    }
}
