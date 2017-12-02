package org.jax.diachromatic.command;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.FragmentFactory;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;

import java.util.ArrayList;
import java.util.List;

import static org.jax.diachromatic.digest.RestrictionEnzyme.parseRestrictionEnzymes;


/**
 * This class coordinates the creation of a digest file that contains a list of restriction fragments across an
 * entire genome.
 */
public class DigestCommand extends Command {
    private static final Logger logger = LogManager.getLogger();
    /** Name of the directory with genome FASTA files and FASTA index (fai) files. */
    private final String genomeDirectoryPath;
    /** Name of output file with list of restriction fragments.*/
    private final String outfilename;
    /** Restriction enzyme that will be used to digest the genome. */
    private RestrictionEnzyme enzyme=null;


    public DigestCommand(String genomeDir,String enzymeName, String outputFile) {
        List<RestrictionEnzyme>  enzymelist = parseRestrictionEnzymes();
        this.enzyme=enzymelist.stream().filter(r->r.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (enzyme==null) {
            logger.fatal(String.format("Could not find restriction enzyme \"%s\". Please correct this and try again",enzymeName ));
            System.exit(1); // todo exception
        }
        genomeDirectoryPath=genomeDir;
        outfilename=outputFile;
    }






    public void execute() {
        FragmentFactory factory=new FragmentFactory(genomeDirectoryPath,outfilename);
        try {
            List<String> enzymes=new ArrayList<>();
            enzymes.add(this.enzyme.getName());
            factory.digestGenome(enzymes);
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }
    }

    @Override
    public String toString() {return "diachromatic:digest";}
}
