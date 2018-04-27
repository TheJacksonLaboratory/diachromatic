package org.jax.diachromatic.command;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.FragmentFactory;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static org.jax.diachromatic.digest.RestrictionEnzyme.parseRestrictionEnzymes;


/**
 * This class coordinates the creation of a digest file that contains a list of restriction fragments across an
 * entire genome.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-19)
 */
public class DigestCommand extends Command {
    private static final Logger logger = LogManager.getLogger();
    /** Name of the directory with genome FASTA files and FASTA index (fai) files. */
    private final String genomeFastaFilePath;
    /** Name of output file with list of restriction fragments.*/
    private final String outfilename;
    /** Restriction enzyme that will be used to digest the genome. */
    private final RestrictionEnzyme enzyme;
    /** size of margin of fragments used for calculating GC and repeat content. */
    private final int marginSize;

    /**
     *
     * @param genomeFile path to the combined FASTA file with all (or all canonical) chromosomes.
     * @param enzymeName name of the enzyme used for digestion (e.g., DpnII)
     * @param outputFile name of output file
     * @param msize margin size (which is used to calculate GC and repeat content)
     */
    public DigestCommand(String genomeFile,String enzymeName, String outputFile, int msize) throws DiachromaticException {
        List<RestrictionEnzyme>  enzymelist = parseRestrictionEnzymes();
        this.enzyme=enzymelist.stream().filter(r->r.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (enzyme==null) {
            throw new DiachromaticException(String.format("Could not find restriction enzyme \"%s\". Please correct this and try again",enzymeName ));
        }
        genomeFastaFilePath =genomeFile;
        if (outputFile.equalsIgnoreCase("default")) {
            String genome=(new File(genomeFile)).getName();
            int i = genome.indexOf("."); // this is nice if the FASTA file has an appropriate name e.g. hg19.fa, but not if the name is e.g. genome.fa
            if (i>0) genome=genome.substring(0,i);
            outfilename=String.format("%s_%s_digest.tsv",genome,enzymeName );
        } else {
            outfilename = outputFile;
        }
        this.marginSize=msize;
    }

    public void execute() {
        FragmentFactory factory=new FragmentFactory(genomeFastaFilePath,outfilename,marginSize);
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
