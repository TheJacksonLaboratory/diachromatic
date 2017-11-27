package org.jax.diachromatic.map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.*;

public class Bowtie2Runner {
    private static final Logger logger = LogManager.getLogger();

    private String pathToBowtie2=null;

    String stdin=null;

    String stderr=null;


    public Bowtie2Runner(String bowtiepath) throws DiachromaticException {
        if (! checkBowtie2(bowtiepath) ){
            throw new DiachromaticException("Could not start bowtie");
        }
    }



    private boolean checkBowtie2(String pathToBowtie2) {
        File bowtie2 = new File(pathToBowtie2);
        if (! bowtie2.exists()) {
            logger.error(String.format("Could not find bowtie2 at %s",pathToBowtie2));
            return false;
        } else if ( ! bowtie2.canExecute()) {
            logger.error(String.format("bowtie2 file at %s was not exectuable",pathToBowtie2));
            return false;
        } else {
            this.pathToBowtie2=bowtie2.getAbsolutePath();
            return true;
        }
    }


    public void run() throws DiachromaticException {
        String args[]=new String[2];
        args[0]=pathToBowtie2;
        args[1]="--help";
        // TODO get rest of command for bowtie2

        String[] dummy=new String[0];
        try {
            // we need to provide the Runtime with the directory in which to start
            // Since we are providing the absolute path of bowtie, we will start in the
            // current working directory (".").
            Process process = Runtime.getRuntime().exec(args,dummy,new File("."));
            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(process.getInputStream()));
            BufferedReader stdError = new BufferedReader(new
                    InputStreamReader(process.getErrorStream()));
            StringBuilder sb = new StringBuilder();
            String s=null;
            while ((s = stdInput.readLine()) != null) {
                sb.append(s+"\n");
            }
            this.stdin=sb.toString();
            sb = new StringBuilder();
            while ((s = stdError.readLine()) != null) {
                sb.append(s+"\n");
            }
            if (sb.length()>0)
                this.stderr=sb.toString();
        } catch (IOException e) {
            String msg = String.format("Could not run bowtie [%s]",e.getMessage());
            throw new DiachromaticException(msg);
        }

        System.out.println("STDOUT="+stdin);
        System.out.println("STDERR="+stderr);
    }



}
