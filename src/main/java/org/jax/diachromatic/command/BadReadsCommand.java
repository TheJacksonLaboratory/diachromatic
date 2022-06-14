package org.jax.diachromatic.command;


import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.BadReadsCounter;
import org.jax.diachromatic.exception.DiachromaticException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;

import java.io.File;
import java.util.concurrent.Callable;

/**
 * Class to coordinate bowtie2-mapping of the truncated FASTQ files followed by Q/C and filtering of the mapped reads.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */
@CommandLine.Command(name = "bad",
        aliases = {"B"},
        mixinStandardHelpOptions = true,
        description = "TODO analyze bad reads.")
public class BadReadsCommand extends Command implements Callable<Integer> {
    private static final Logger logger = LoggerFactory.getLogger(BadReadsCommand.class);

    @CommandLine.Option(names={"-r", "--rejected-pairs-bam"}, required = true, description = "Path to BAM file with unique valid pairs produced using the align command.", order = 3)
    private String rejectedPairsBamFile;

    /** Path to the genome digest file produced by GOPHER. */
    @CommandLine.Option(names={"-d","--digest-file"}, required = true, description = "Path to GOPHER digest file.", order = 4)
    private String digestFile;




    @Override
    public Integer call() throws DiachromaticException {
       // input BAM file with HTS-JDK
        // iterate over (rejected) read-pairs
       // count up read pairs per bait, whereby we have 5D -- B -- 3D, where 5D is the digest
        // 5' to the bait (B) and 3D is the 3' digect
        //

        makeOutdirectoryIfNeeded();

        logger.trace(String.format("About to read rejected read pairs file from %s", rejectedPairsBamFile));
        DigestMap digestMap = new DigestMap(digestFile);
        String outputName=String.format("%s%s%s.tsv", outputDir, File.separator,filenamePrefix);
        SamReader reader = SamReaderFactory.makeDefault().open(new File(rejectedPairsBamFile));
        BadReadsCounter counter = new BadReadsCounter(reader, digestMap);
        counter.countInteractions();

        counter.writeToFile(outputName);

        // todo -- BAM tags., e.g. 	YY:Z:UL

        // Un-ligated due to size (Tag: UL) ZAEHLT, die anderen nicht
        // Too short chimeric (Tag: TS) -- nachlesen, wohl ja

        return 0;
    }

}
