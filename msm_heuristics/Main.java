import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import MSMMean.MsmMean;
import MSMMean.Pair;
import ThreeMsmMedian.ThreeMsmMedian;
import buckets_median.AverageBucketsMsmMedian;
import buckets_median.BucketsMsmMedian;
import buckets_median.RandomBucketsMergeMsmMedian;
import buckets_median.RandomBucketsMsmMedian;
import buckets_median.RandomMsmMedian;
import evolutionary_algorithm.BucketMsmMedianEA;
import evolutionary_algorithm.GenericMsmMedianEA;

public class Main {
    
    static DecimalFormat df = new DecimalFormat("#.######");

    static int k = 4;
    static int n = 10;
    static int amount = 3;

    public static void main(String[] args) {
        int[] bucketsize = {2, 3, 4, 5};
        int[] bucketradius = {1, 2, 3};
        int[] amountRandom = {10_000, 100_000, 1_000_000};
        int iteration = 10;
        int[] generations = {10_000, 100_000, 1_000_000};

        String currentLine;
        BufferedReader reader;
        
        double[][] timeseries;
        
        try {
            String[] sets = {
                "50words", 
                "Adiac", 
                "Beef", 
                "CBF", 
                "Coffee", 
                "ECG200", 
                "FaceAll", 
                "FaceFour", 
                "FISH", 
                "Gun_Point",
                "Lighting2", 
                "Lighting7", 
                "OliveOil", 
                "OSULeaf", 
                "SwedishLeaf", 
                "synthetic_control", 
                "Trace", 
                "Two_Patterns", 
                "wafer", 
                "yoga"
            };

            System.out.println();
            while(k <= 4) {
                while(n <= 17) {

                    File csvFile = new File("experiments/randomMsmMedian/heuristic_k_" + k + "_n_" + n +"_only_random_msm.csv");
                    PrintWriter pw = new PrintWriter(csvFile);
                    pw.println("heuristic;k;n;dataset;c;mean;cost;time in ms;");
                    // header

                    for(int currentSetIndex = 0; currentSetIndex < sets.length; currentSetIndex++) {
                        String file = "/Users/hieu/Library/CloudStorage/OneDrive-Personal/Studium/BA Heuristische Verfahren für die Berechnung von MSM Median/randomSamplesUCRArchive/k_" + k + "/n_"+ n 
                        +"/mixedClasses/" + sets[currentSetIndex] + "_TRAIN_k_"+ k + "_n_" + n + ".csv";
                        System.out.println("Computing for k = " + k + ", n = " + n + ": ");
                        System.out.println("with the current file: ");
                        System.out.println(file);

                        double c = getC(file);
                        if(c == -1) {
                            throw new IllegalArgumentException("value of c is wrong");
                        }

                        reader = new BufferedReader(new FileReader(file));
                        int index = 0;
                        timeseries = new double[k][];
                        while((currentLine = reader.readLine()) != null) {
                            String[] sArr = currentLine.split(",");
                            double[] ts = Arrays.stream(sArr)
                                    .mapToDouble(Double::parseDouble)
                                    .toArray();
                            timeseries[index] = ts;
                            index++;
                        }

                        // for(double[] ts: timeseries) {
                        //     System.out.println(Arrays.toString(ts));
                        // }

                        // msm mean exact
                        for(int i = 0; i < amount; i++) {
                            // msmMeanExact(timeseries, c, sets[currentSetIndex],pw);
                        }
                        // 3 msm mean
                        for(int i = 0; i < amount; i++) {
                            // threeMsmMean(timeseries, c, sets[currentSetIndex], pw);
                        }
                        // bucket msm mean
                        for(int j = 0; j < bucketsize.length; j++) {
                            for(int i = 0; i < amount; i++) {
                                // bucketsMsmMean(timeseries, c, bucketsize[j], sets[currentSetIndex], pw);
                            }
                        } 
                        // random msm median
                        for(int r = 0; r < amountRandom.length; r++) {
                            for(int i = 0; i < amount; i++) {
                                randomMsmMedian(timeseries, c, amountRandom[r], sets[currentSetIndex], pw);
                            }
                        }
                        // average bucket msm mean
                        for(int j = 0; j < bucketsize.length; j++) {
                            for(int r = 0; r < amountRandom.length; r++) {
                                for(int i = 0; i < amount; i++) {
                                    // averageBucketsMsmMean(timeseries, c, amountRandom[r], bucketsize[j], iteration, sets[currentSetIndex], pw);
                                }
                            }
                        }
                        
                        // random buckets msm mean
                        for(int j = 0; j < bucketsize.length; j++) {
                            for(int r = 0; r < amountRandom.length; r++) {
                                for(int i = 0; i < amount; i++) {
                                    // randomBucketsMsmMean(timeseries, c, amountRandom[r], bucketsize[j], sets[currentSetIndex], pw);
                                }
                            }
                        }
                        
                        // random buckets merge msm mean
                        for(int j = 0; j < bucketsize.length; j++) {
                            for(int r = 0; r < amountRandom.length; r++) {
                                for(int i = 0; i < amount; i++) {
                                    // randomBucketsMergeMsmMean(timeseries, c, amountRandom[r], bucketsize[j], 2*c, sets[currentSetIndex], pw);
                                }
                            }
                        }
                        // generic msm mean ea
                        for(int g = 0; g < generations.length; g++) {
                            for(int i = 0; i < amount; i++) {
                                // genericMsmMeanEA(timeseries, c, generations[g], sets[currentSetIndex], pw);
                            }
                        }

                        // buckets msm mean ea
                        for(int j = 0; j < bucketradius.length; j++) {
                            for(int g = 0; g < generations.length; g++) {
                                for(int i = 0; i < amount; i++) {
                                    // bucketMsmMeanEA(timeseries, c, bucketradius[j], generations[g], sets[currentSetIndex], pw);
                                }
                            }
                        }

                        

                        
                    }
                    pw.close();

                    n++;
                }
                n = 10;
                k++;
            }
            
        } 
        catch (Exception e) {
            System.out.println();
            System.out.println("No file or file not found.");
        }

    }

    private static double getC(String file) {
        String[] s1 = {"50words", "Adiac", "ECG200", "FaceAll", "FaceFour", "Lighting7", "SwedishLeaf", "Two_Patterns", "wafer"};
        String[] s2 = {"Beef", "CBF", "FISH", "OSULeaf", "synthetic_control", "yoga"};
        String[] s3 = {"Coffee", "Gun_Point", "Lighting2", "OliveOil", "Trace"};
        if(Arrays.stream(s1).anyMatch(file::contains)) {
            return 1.;
        } 
        else if(Arrays.stream(s2).anyMatch(file::contains)){
            return 0.1;
        }
        else if(Arrays.stream(s3).anyMatch(file::contains)) {
            return 0.01;
        } 
        else {
            return -1.;
        }
    }

    private static void line() {
        String line = "---------------------------------";
        System.out.println();
        System.out.println(line);
    }



    private static void msmMeanExact(double[][] timeseries, double c, String currentSet,PrintWriter pw) {
        String methodName = "MSM-Mean Exact";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        MsmMean msmMean = new MsmMean(timeseries, c);
        Pair<double[], Double> msm = msmMean.calcMean();
        System.out.println(Arrays.toString(msm.getFirst()));
        System.out.println("Cost: " + msm.getSecond());
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(msm.getFirst())+";"+df.format(msm.getSecond())+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }
    
    private static void threeMsmMean(double[][] timeseries, double c, String currentSet, PrintWriter pw) {    
        String methodName = "3-MSM-Mean";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        ThreeMsmMedian bmsmMean = new ThreeMsmMedian(timeseries, c);
        double[][] result = bmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void bucketsMsmMean(double[][] timeseries, double c, int bucketSize, String currentSet, PrintWriter pw) {
        String methodName = "Buckets MSM-Mean";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        BucketsMsmMedian bmsmMean = new BucketsMsmMedian(timeseries, c, bucketSize);
        double[][] result = bmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": bucketsize:"+bucketSize+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void averageBucketsMsmMean(double[][] timeseries, double c, int amountRandom, int bucketSize, int iteration, String currentSet, PrintWriter pw) {
        String methodName = "Average Buckets MSM-Mean";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        AverageBucketsMsmMedian bmsmMean = new AverageBucketsMsmMedian(timeseries, c, amountRandom, bucketSize, iteration);
        double[][] result = bmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": bucketsize:"+bucketSize+", rand:"+amountRandom+", iteration:"+iteration+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void randomMsmMedian(double[][] timeseries, double c, int amountRandom, String currentSet, PrintWriter pw) {
        String methodName = "Random MSM Median";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        RandomMsmMedian rmsmMean = new RandomMsmMedian(timeseries, c, amountRandom);
        double[][] result = rmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": rand:"+amountRandom+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void randomBucketsMsmMean(double[][] timeseries, double c, int amountRandom, int bucketSize, String currentSet, PrintWriter pw) {
        String methodName = "Random Buckets MSM-Mean";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        RandomBucketsMsmMedian bmsmMean = new RandomBucketsMsmMedian(timeseries, c, amountRandom, bucketSize);
        double[][] result = bmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": bucketsize:"+bucketSize+", rand:"+amountRandom+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void randomBucketsMergeMsmMean(double[][] timeseries, double c, int amountRandom, int bucketSize, double mergeRange, String currentSet, PrintWriter pw) {
        String methodName = "Random Buckets Merge MSM-Mean";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        RandomBucketsMergeMsmMedian bmsmMean = new RandomBucketsMergeMsmMedian(timeseries, c, amountRandom, bucketSize, mergeRange);
        double[][] result = bmsmMean.calc();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": bucketsize:"+bucketSize+", rand:"+amountRandom+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void genericMsmMeanEA(double[][] timeseries, double c, int generations, String currentSet, PrintWriter pw) {
        String methodName = "Generic MSM-Mean EA";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        GenericMsmMedianEA bmsmMean = new GenericMsmMedianEA(timeseries, c, generations);
        double[][] result = bmsmMean.evolve();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": generations:"+generations+";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }

    private static void bucketMsmMeanEA(double[][] timeseries, double c, int bucketRadius,int generations, String currentSet, PrintWriter pw) {
        String methodName = "Buckets MSM-Mean EA";
        long startTime = System.nanoTime();
        System.out.println(methodName);
        BucketMsmMedianEA bmsmMean = new BucketMsmMedianEA(timeseries, c, bucketRadius,generations);
        double[][] result = bmsmMean.evolve();
        System.out.println(Arrays.toString(result[0]));
        System.out.println("Cost: " + result[1][0]);
        long endTime = System.nanoTime();
        long totalTime = endTime - startTime;
        System.out.println("Running Time: " + TimeUnit.NANOSECONDS.toMillis(totalTime) + " ms");
        line();
        pw.println(methodName+": generations:"+generations+", bucketradius: "+bucketRadius +";"+k+";"+n+";"+currentSet+";"+c+";"+Arrays.toString(result[0])+";"+df.format(result[1][0])+";"+TimeUnit.NANOSECONDS.toMillis(totalTime));
    }
}
