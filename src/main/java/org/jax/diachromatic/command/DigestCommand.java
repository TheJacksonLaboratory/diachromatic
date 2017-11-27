package org.jax.diachromatic.command;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.FragmentFactory;
import org.jax.diachromatic.exception.DiachromaticException;

import java.util.ArrayList;
import java.util.List;

public class DigestCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    private String genomeDirectoryPath=null;
    /** Name of output file. TODO set this dynamically */
    private String outfilename="hicupCloneDigest.txt";


    public DigestCommand(String genomeDir) {
        genomeDirectoryPath=genomeDir;
    }

    public DigestCommand(String genomeDir,String outputFile) {

    }


    public void execute() {
        FragmentFactory factory=new FragmentFactory(genomeDirectoryPath,outfilename);
        try {
            factory.indexFASTAfilesIfNeeded();
            String testEnzyme = "DpnII";
            List<String> enzymes=new ArrayList<>();
            enzymes.add(testEnzyme);
            factory.digestGenome(enzymes);
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }
        // implement digesting a genome
    }

    @Override
    public String toString() {return "diachromatic:digest";}
}
