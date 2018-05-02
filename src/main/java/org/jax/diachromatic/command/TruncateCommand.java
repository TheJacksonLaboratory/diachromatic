package org.jax.diachromatic.command;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.Diachromatic;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.truncation.Truncator;

import java.util.Collections;
import java.util.List;

import static org.jax.diachromatic.digest.RestrictionEnzyme.parseRestrictionEnzymes;

/**
 * This command parallels the HiCUP truncate Perl script with various modifcations. The following description is
 * based on the description from the Babraham bioinformatics HiCUP documentation.
 * <P>
 *     Valid Hi-C pairs are chimeric reads that are made up of fragments from two different regions of the genome.
 *     With capture Hi-C, a typical valid read pair might comprise a DNA sequence from a promoter and a DNA sequence
 *     from an enhancer that is regulating the promoter. In most cases, the one of the reads of the read paier will map
 *     to a single ligation fragment, and the reverse read will map to another fragment. However, this is not always true
 *     because the Hi-C ligation junction cvan be located within one of the sequenced reads. The truncater attempts to
 *     address this situation (which could lead to the read with the Hi-C junction not being mapped during the mapping
 *     step), buy deleting sequenced that is downstream of the enzyme recognition site. For example, if the forward read
 *     is entirely contained withint one ligation fragment and the reverse read starts in another fragment, leads into the
 *     ligation junction, and then continues and finishes with part of the fragment of the forward read, then the
 *     truncation step will remove the part of the reverse read that maps to the first ligation fragment.
 * </P>
 * <P>
 *     The method requires the names of the files to be processed as well as the restriction enzymes used to process the
 *     DNA. We will expect that there are two FAASTQ files that are gzipped that represent the paired end sequences.
 * </P>
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.2 (2018-01-05)
 */
public class TruncateCommand extends Command {
    private static final Logger logger = LogManager.getLogger();
    private final String outdir;

    private final String fastaqFile1;
    private final String fastaqFile2;

    private String truncatedFastaqFile1=null;
    private String truncatedFastaqFile2=null;

    private Truncator truncator = null;

    private RestrictionEnzyme re=null;

    public TruncateCommand (String file1, String file2, String enzymeName, String outdir, String outprefix) throws DiachromaticException {
        System.out.println(outdir);
        System.out.println(outprefix);
        this.outdir=outdir;
        this.fastaqFile1=file1;
        this.fastaqFile2=file2;
        List<RestrictionEnzyme>  enzymelist = parseRestrictionEnzymes();
        re=enzymelist.stream().filter(r->r.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (re==null) {
            throw new DiachromaticException(String.format("Could not identify restriction enzyme for \"%s\"",enzymeName));
        }
        truncator = new Truncator(fastaqFile1,fastaqFile2,re,outdir,outprefix);
    }

    public void execute() {
        logger.trace(String.format("Starting truncate command on files %s and %s",fastaqFile1,fastaqFile2 ));
        try {
            truncator.parseFASTQ();
        } catch (DiachromaticException e) {
            logger.fatal("Error encountered while truncating FASTQ reads: ",e);
        }
    }

    @Override
    public String toString() {return "diachromatic:truncate";}
}
