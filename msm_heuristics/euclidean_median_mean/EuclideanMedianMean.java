package euclidean_median_mean;


import java.util.Arrays;
import java.util.concurrent.TimeUnit;

/**
 * Creates a time series which represents the mean of given time series calculated with the 
 * euclidean distance. 
 * Calculates the cost between mean and all time series. There are only move costs at the beginning.
 * At the end, if the lengths of the time series are not equal, there are side merges with "small" moves.
 */
public class EuclideanMedianMean {

    protected final double[][] timeseries;
    protected final int k; 
    protected final double c;
    protected double[] medianMean;

    protected int minDimension = Integer.MAX_VALUE;

    public EuclideanMedianMean(double[][] timeseries, double c) {
        this.timeseries = timeseries;
        this.k = timeseries.length;
        this.c = c;

        for(double[] ts : this.timeseries) {
            this.minDimension = Math.min(minDimension, ts.length);
        }
        this.medianMean = new double[minDimension];
        
        for(int i = 0; i < minDimension; i++) {
            double[] currentColumn = new double[k];
            for(int j = 0; j < k; j++) {
                currentColumn[j] = timeseries[j][i];
            }
            medianMean[i] = getMedian(currentColumn);
        }
            
    }

    /**
     * Returns you a median of type double. We differentiate between two cases: the length of the 
     * given array can be even or odd.
     * @param arr given array
     * @return returns the median of type double
     */
    public double getMedian(double[] arr) {
        Arrays.sort(arr);
        double median;
        if (arr.length % 2 == 0)
            median = ((double)arr[arr.length/2] + (double)arr[arr.length/2 - 1])/2;
        else
            median = (double) arr[arr.length/2];
        return median;
    }

    /**
     * Returns the cost of the move operation.
     * @param x an element x of a time series 
     * @param y an element y of a time series
     * @return the cost of type double
     */
    private double costMove(double x, double y) {
        return Math.abs(y - x);
    }

    /**
     * Returns the cost of the side merges and the small moves.
     * @param lastMedianMean gives the last element of the mean.
     * @return the cost of type doble
     */
    private double costSideMerge(double lastMedianMean) {
        double cost = 0.;
        for(int i = 0; i < this.k; i++) {
            for(int j = this.medianMean.length; j < this.timeseries[i].length;j++) {
                cost += this.c + Math.abs(this.timeseries[i][j] - this.timeseries[i][j-1]);
            }
        }
        return cost;
    }

    /**
     * Calculates the cost between the mean and all time series.
     */
    private void calcCost() {
        double sum = 0;
        for(int i = 0; i < this.medianMean.length; i++) {
            for(int j = 0; j < this.k; j++) {
                sum += costMove(this.medianMean[i], this.timeseries[j][i]);
                // System.out.println("Cost for current position of Median (" + (i+1) + ") and current position of time series (" + (j+1) + ", " + (i+1) + "): " + costMove(this.medianMean[i], this.timeseries[j][i]));
            }
        }

        sum  = sum + costSideMerge(this.medianMean[medianMean.length-1]);
        System.out.println("Total cost: " + sum);
    }



    public static void main(String[] args) {
        double c = 0.1;

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
        System.out.println("Computing for k = " + k + " and c = " + c);
        System.out.println();

        System.out.println("Median Mean: ");
        long startTime = System.nanoTime();
        EuclideanMedianMean euclideanMedianMean = new EuclideanMedianMean(timeseries, c);
        System.out.println(Arrays.toString(euclideanMedianMean.medianMean));
        euclideanMedianMean.calcCost();
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
    }
}