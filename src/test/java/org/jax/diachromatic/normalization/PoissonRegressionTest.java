package org.jax.diachromatic.normalization;

//import com.github.chen0040.data.frame.BasicDataFrame;
//import com.github.chen0040.data.frame.DataRow;
//import com.github.chen0040.glm.enums.GlmDistributionFamily;
//import com.github.chen0040.glm.solvers.Glm;
//import com.github.chen0040.glm.enums.GlmSolverType;
//import com.github.chen0040.data.frame.DataFrame;
//import org.apache.logging.log4j.LogManager;
//import org.jax.diachromatic.normalization.PoissonRegression;
//import org.apache.logging.log4j.Logger;


import org.junit.Test;

import java.io.*;

public class PoissonRegressionTest {
    //private static final Logger logger = LogManager.getLogger();



    /*
       The following code is intended to reproduce an example for Poisson Regression Model for Count Data in R
       from here: https://onlinecourses.science.psu.edu/stat504/node/169

       Using a java GLM package from here: https://github.com/chen0040/java-glm

       See also issue #21
     */

    @Test
    public void reproduceHorseshoeCrabExample() {

        /* Read data provided for the tutorial: The structure of the file crab.txt provided by the tutorial is as follows:

               crab=read.table("crab.txt")
               colnames(crab)=c("Obs","C","S","W","Wt","Sa")

        */

        //DataFrame dataFrameCrab = new BasicDataFrame();

        BufferedReader br = null;
        try {

            br = new BufferedReader(new FileReader("src/test/resources/data/glm_test_data/crab.txt"));
            String line;

            while ((line = br.readLine()) != null) {

                String line2=line.trim();
                String[] A = line2.split(" +");

                Integer observation=Integer.parseInt(A[0]);
                String color=A[1];
                Integer spine_condition=Integer.parseInt(A[2]);
                Double weight=Double.parseDouble(A[3]);
                Double width=Double.parseDouble(A[4]);
                Integer satelite_number=Integer.parseInt(A[5]);
                //System.out.println(observation + "\t " + color + "\t " + spine_condition + "\t " + weight + "\t " + width + "\t " + satelite_number);

                //DataRow row = dataFrameCrab.newRow();
                //row.setCell("observation", observation);
                //row.setCategoricalCell("color", color);
                //row.setCell("spineCondition", spine_condition);
                //row.setCell("weight", weight);
                //row.setCell("width", width);
                //row.setTargetCell("sateliteNumber", satelite_number);
                //dataFrameCrab.addRow(row);

            }
            //dataFrameCrab.lock();

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) {
                    br.close();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }

        /* Fit a Poisson regression model with only one predictor: width (W) */

        System.out.println("Set up model");
/*
        Glm glm = Glm.linear();
        glm.setDistributionFamily(GlmDistributionFamily.Poisson);
        glm.setSolverType(GlmSolverType.GlmIrls);

        System.out.println("Distribution family: " + glm.getDistributionFamily());
        System.out.println("All columns: " + dataFrameCrab.getAllColumns());
        System.out.println("Input columns: " + dataFrameCrab.getInputColumns());
        System.out.println("Output column: " + dataFrameCrab.getOutputColumns().get(0).getColumnName());
        System.out.println("Levels: " + dataFrameCrab.getLevels());

        glm.fit(dataFrameCrab);

        System.out.println(glm.getCoefficients());
*/
    }




}