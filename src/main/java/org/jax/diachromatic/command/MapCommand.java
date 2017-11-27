package org.jax.diachromatic.command;

import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.map.Bowtie2Runner;
import org.jax.diachromatic.map.SAMPairer;

import java.io.IOException;

public class MapCommand extends Command {

    private final String bowtiepath;


    private String pathToBowtieIndex=null;

    private String pathToInputFastq1 =null;

    private String pathToInputFastq2 =null;

    private String outname=null;

    public MapCommand(String bowtie, String btIndexPath, String inputFastqPath1,String inputFastqPath2, String outnam) {
        bowtiepath=bowtie;
        pathToBowtieIndex=btIndexPath;
        pathToInputFastq1 =inputFastqPath1;
        pathToInputFastq2 =inputFastqPath2;
        outname=outnam;

    }

    public void execute() {
        String outname1="tempoutname1.sam";
        String outname2="tempoutname2.sam";

        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath,pathToBowtieIndex, pathToInputFastq1,outname1);
            runner.run();
            Bowtie2Runner runner2 = new Bowtie2Runner(bowtiepath,pathToBowtieIndex, pathToInputFastq2,outname2);
            runner2.run();
            SAMPairer pairer = new SAMPairer(outname1,outname2);
            pairer.pair();
        } catch (DiachromaticException e){
            e.printStackTrace();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }
    }
    @Override
    public String toString() {return "diachromatic:map";}
}
