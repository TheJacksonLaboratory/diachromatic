package org.jax.diachromatic.command;


import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.Counter;
import org.jax.diachromatic.exception.DiachromaticException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.concurrent.Callable;

/**
 * Class to coordinate counting of valid read pairs between pairs of restriction digests.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */

@CommandLine.Command(name = "count",
        aliases = {"C"},
        mixinStandardHelpOptions = true,
        description = "count  valid pairs between pairs of restriction fragments from a BAM file creted in the alignstep with a GOPHER digest file.")
public class CountCommand extends Command implements Callable<Integer> {
    private static final Logger logger = LoggerFactory.getLogger(CountCommand.class);

    //@CommandLine.Option(names = {"-j", "--jannovar"}, description = "prefix for output files (default: ${DEFAULT-VALUE} )")
    /** Path to BAM file containing unique valid pairs. */
    @CommandLine.Option(names={"-v", "--valid-pairs-bam"}, required = true, description = "Path to BAM file with unique valid pairs produced using the align command.", order = 3)
    private String validPairsBamFile = null;

    /** Path to the genome digest file produced by GOPHER. */
    @CommandLine.Option(names={"-d","--digest-file"}, required = true, description = "Path to GOPHER digest file.", order = 4)
    private String digestFile = null;

    /** Aggregate counts for different read pair orienations. */
    @CommandLine.Option(names={"-s", "--split-counts"},description = "Split counts for different read pair orientations.", order = 5)
    private boolean split=false;


    public CountCommand() {
    }

    @Override
    public Integer call() throws DiachromaticException {

        makeOutdirectoryIfNeeded();

        logger.trace(String.format("About to read digests from %s",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);

        String outputDirAndFilePrefix=String.format("%s%s%s", outputDir, File.separator,filenamePrefix);

        SamReader reader = SamReaderFactory.makeDefault().open(new File(validPairsBamFile));

        Counter counter = new Counter(reader, digestMap, outputDirAndFilePrefix, split);
        try {
            logger.trace("About to determine interaction counts...");
            counter.countInteractions();
            logger.trace("...done with counting!");
            logger.trace("About to print the results...");
            counter.printInteractionCountsMapAsCountTable();
            counter.printInteractionCountsMapInWashUSimpleTextFormat();
            counter.printFragmentInteractionCountsMapAsCountTable();
            counter.printStatistics();
            logger.trace("...done!");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        return 0;
    }
    @Override
    public String toString() {return "diachromatic:count";} //???
}
