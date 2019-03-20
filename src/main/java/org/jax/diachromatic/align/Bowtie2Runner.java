package org.jax.diachromatic.align;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.*;
import java.util.Arrays;
import java.util.stream.Collectors;

public class Bowtie2Runner {
    private static final Logger logger = LogManager.getLogger();

    private String pathToBowtie2;

    private final String pathToBowtieIndex;

    private final String pathToInputFastq;

    private final String outname;

    private final int threadNum;

    private String stdin = null;

    private String stderr = null;


    public Bowtie2Runner(String bowtiepath, String btIndexPath, String inputFastqPath, String outnam, Integer threadNum) throws DiachromaticException {
        if (! checkBowtie2(bowtiepath) ){
            throw new DiachromaticException("Could not start bowtie");
        }

        pathToBowtieIndex=btIndexPath;
        pathToInputFastq=inputFastqPath;
        outname=outnam;
        checkExistenceOfInputFile();
        this.threadNum=threadNum;
    }


    /** Throw an exception of the input file (which should be something like
     * foo.truncated_R1.fastq.gz) cannot be found
     * @throws DiachromaticException if the input file cannot be found.
     */
    private void checkExistenceOfInputFile() throws  DiachromaticException{
        File f = new File(pathToInputFastq);
        if (!f.exists()) {
            throw new DiachromaticException("Could not find input truncated FASTQ file at " + pathToInputFastq);
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

    /**
     * Run bowtie.
     *
     * @throws DiachromaticException
     */
    public void run() throws DiachromaticException {
        String[] args = new String[11];
        args[0]=pathToBowtie2;
        args[1]="--very-sensitive";
        //args[2]="--no-unal";
        args[2]="-p";
        args[3]=String.valueOf(threadNum);
        args[4]="--reorder"; // keep same order of records as in FASTQ
        args[5]="-x";
        args[6]=pathToBowtieIndex;
        args[7]="-U"; // unpaired reads to be aligned
        args[8]=pathToInputFastq; // Input FASTQ file (just one!)
        args[9]="-S";
        args[10]=outname;// summarize name


        String btcomd= Arrays.stream(args).collect(Collectors.joining(" "));
        logger.trace("Running: "+btcomd);

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
