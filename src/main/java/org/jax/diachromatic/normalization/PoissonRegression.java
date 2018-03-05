package org.jax.diachromatic.normalization;

import java.util.Arrays;

import org.apache.spark.ml.regression.GeneralizedLinearRegression;
import org.apache.spark.ml.regression.GeneralizedLinearRegressionModel;
import org.apache.spark.ml.regression.GeneralizedLinearRegressionTrainingSummary;
import org.apache.spark.sql.Dataset;
import org.apache.spark.sql.Row;
import org.apache.spark.sql.SparkSession;
import org.apache.spark.sql.types.DataTypes;
import org.apache.spark.sql.types.Metadata;
import org.apache.spark.sql.types.StructField;
import org.apache.spark.sql.types.StructType;


public class PoissonRegression {

    private final SparkSession spark;
    public PoissonRegression() {
        spark = SparkSession
                .builder()
                .appName("Diachromatic")
               // .config("spark.some.config.option", "some-value")
                .getOrCreate();
    }

    private void init() {


    }
}
