package org.jax.diachromatic.command;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.FragmentFactory;

public class DigestCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    private String genomeDirectoryPath=null;


    public DigestCommand(String genomeDir) {
        genomeDirectoryPath=genomeDir;
    }

    public DigestCommand(String genomeDir,String outputFile) {

    }


    public void execute() {
        FragmentFactory factory=new FragmentFactory(genomeDirectoryPath);
        factory.indexFASTAfilesIfNeeded();
        String testEnzyme = "DpnII";
        factory.cutWithEnzyme(testEnzyme);
        // implement digesting a genome
    }
}
