package org.jax.diachromatic.command;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.FragmentFactory;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;

import java.util.ArrayList;
import java.util.List;

import static org.jax.diachromatic.digest.RestrictionEnzyme.parseRestrictionEnzymes;

public class DigestCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    private String genomeDirectoryPath=null;
    /** Name of output file. TODO set this dynamically */
    private String outfilename="hicupCloneDigest.txt";

    private RestrictionEnzyme enzyme=null;


    public DigestCommand(String genomeDir,String restrictionEnzyme) {
        init(restrictionEnzyme);
        genomeDirectoryPath=genomeDir;
    }

    public DigestCommand(String genomeDir,String restrictionEnzyme, String outputFile) {
        init(restrictionEnzyme);
        genomeDirectoryPath=genomeDir;
        outfilename=outputFile;
    }


    private void init(String enzymeName) {
        List<RestrictionEnzyme>  enzymelist = parseRestrictionEnzymes();
        this.enzyme=enzymelist.stream().filter(r->r.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (enzyme==null) {
            logger.fatal(String.format("Could not find restriction enzyme \"%s\". Please correct this and try again",enzymeName ));
            System.exit(1); // todo exception
        }
    }




    public void execute() {
        FragmentFactory factory=new FragmentFactory(genomeDirectoryPath,outfilename);
        try {
            factory.indexFASTAfilesIfNeeded();
            List<String> enzymes=new ArrayList<>();
            enzymes.add(this.enzyme.getName());
            factory.digestGenome(enzymes);
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }
        // implement digesting a genome
    }

    @Override
    public String toString() {return "diachromatic:digest";}
}
