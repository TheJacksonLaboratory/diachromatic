package org.jax.diachromatic.command;

import com.google.common.collect.ImmutableList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.map.Bowtie2Runner;
import org.jax.diachromatic.map.Digest;
import org.jax.diachromatic.map.SAMPairer;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class MapCommand extends Command {
    private static final Logger logger = LogManager.getLogger();
    private final String bowtiepath;


    private String pathToBowtieIndex=null;

    private String pathToInputFastq1 =null;

    private String pathToInputFastq2 =null;

    private String outname=null;

    private String digestFile=null;

    public MapCommand(String bowtie, String btIndexPath, String inputFastqPath1,String inputFastqPath2, String outnam, String digest) {
        bowtiepath=bowtie;
        pathToBowtieIndex=btIndexPath;
        pathToInputFastq1 =inputFastqPath1;
        pathToInputFastq2 =inputFastqPath2;
        outname=outnam;
        digestFile=digest;

    }

    public void execute() {
        String outname1="tempoutname1.sam";
        String outname2="tempoutname2.sam";
        logger.trace(String.format("About to read digests from %s",digestFile ));
        Map<String,List<Digest>> digestmap = Digest.readDigests(digestFile);
        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath,pathToBowtieIndex, pathToInputFastq1,outname1);
            runner.run();
            Bowtie2Runner runner2 = new Bowtie2Runner(bowtiepath,pathToBowtieIndex, pathToInputFastq2,outname2);
            runner2.run();
            SAMPairer pairer = new SAMPairer(outname1,outname2,digestmap);
            pairer.pair();
            pairer.printStatistics();
        } catch (DiachromaticException e){
            e.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

    }
    @Override
    public String toString() {return "diachromatic:map";}
}
