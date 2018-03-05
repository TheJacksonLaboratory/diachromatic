package org.jax.diachromatic.normalization;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.spark.ml.linalg.VectorUDT;
import org.apache.spark.ml.linalg.Vectors;
import org.apache.spark.ml.regression.GeneralizedLinearRegression;
import org.apache.spark.ml.regression.GeneralizedLinearRegressionModel;
import org.apache.spark.ml.regression.GeneralizedLinearRegressionTrainingSummary;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.RowFactory;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.types.*;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.map.SAMPairerTest;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * We are testing the code using data from here
 * https://onlinecourses.science.psu.edu/stat504/book/export/html/169
 * For the "width" test, we regress the Satellites against the carapace width (W)
 * See the data in src/test/resources/data
 * The column names are colnames(crab)=c("Obs","C","S","W","Wt","Sa")
 *  female crab had any other males, called satellites, residing near her. Explanatory variables that are thought to
 *  affect this included the female crab’s color (C), spine condition (S), weight (Wt), and carapace width (W).
 *  According to the R implementation, The estimated model is: $log (\hat{\mu_i})$ = -3.30476 + 0.16405Wi
 *  <p>
 *      Currently trying to understand the syntax of Poisson regression in Apache spark
 *  </p>
 */
public class PoissonRegressionTest {

    private static List<Double> satellites;
    private static List<Double> width;

    private static final double EPSILON=0.001;

    @BeforeClass
    public static void init() throws DiachromaticException, IOException  {
        ClassLoader classLoader = SAMPairerTest.class.getClassLoader();
        String horseshoeCrabData = classLoader.getResource("data/crab.txt").getFile();
        satellites=new ArrayList<>();
        width=new ArrayList<>();
        BufferedReader br = new BufferedReader(new FileReader(horseshoeCrabData));
        String line;
        while ((line=br.readLine())!=null) {
            String A[]= line.split("\\s+");
            double n_satellites=Double.parseDouble(A[5]);
            double wdth = Double.parseDouble(A[3]);
            satellites.add(n_satellites);
            width.add(wdth);
        }
    }



    @Test
    public void sparkTest() {
        SparkSession spark = SparkSession
                .builder()
                .appName("Diachromatic")
                .config("spark.master", "local")
                .getOrCreate();
        StructType schema = new StructType(new StructField[]{
                new StructField("satellites",
                        DataTypes.DoubleType, false,
                        Metadata.empty()),
                new StructField("width",
                        new VectorUDT(), false,
                        Metadata.empty())});
        List<Row> rows = new ArrayList<>();
        for (int i=0;i<satellites.size();++i) {
            double s = satellites.get(i);
            double w = width.get(i);
            rows.add(RowFactory.create(s,Vectors.dense(w)));
        }
        System.out.println("ROWS      " +rows.toString());
        System.out.println("SCHEME      " +schema.toString());
        Dataset<Row> df = spark.createDataFrame(rows, schema);
        GeneralizedLinearRegression glr = new GeneralizedLinearRegression()
                .setFamily("poisson")
                .setLink("log")
                .setFeaturesCol("width")
                .setLabelCol("satellites")
                .setMaxIter(1)
                .setRegParam(0.0);

// Fit the model
        GeneralizedLinearRegressionModel model = glr.fit(df);

// Print the coefficients and intercept for generalized linear regression model
        System.out.println("Coefficients: " + model.coefficients());
        System.out.println("Intercept: " + model.intercept());

// Summarize the model over the training set and print out some metrics
        GeneralizedLinearRegressionTrainingSummary summary = model.summary();
        System.out.println("Coefficient Standard Errors: "
                + Arrays.toString(summary.coefficientStandardErrors()));
        System.out.println("T Values: " + Arrays.toString(summary.tValues()));
        System.out.println("P Values: " + Arrays.toString(summary.pValues()));
        System.out.println("Dispersion: " + summary.dispersion());
        System.out.println("Null Deviance: " + summary.nullDeviance());
        System.out.println("Residual Degree Of Freedom Null: " + summary.residualDegreeOfFreedomNull());
        System.out.println("Deviance: " + summary.deviance());
        System.out.println("Residual Degree Of Freedom: " + summary.residualDegreeOfFreedom());
        System.out.println("AIC: " + summary.aic());
        System.out.println("Deviance Residuals: ");
        summary.residuals().show();
    }
}
