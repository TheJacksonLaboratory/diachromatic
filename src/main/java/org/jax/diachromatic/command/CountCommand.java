package org.jax.diachromatic.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.Counter;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;
import java.io.FileNotFoundException;

/**
 * Class to coordinate counting of valid read pairs between pairs of restriction digests.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */
@Parameters(commandDescription = "The count command takes a BAM file containing unique valid read pairs determined during the align step as well as a GOPHER digest file and counts valid pairs between pairs of restriction fragments.")
public class CountCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to BAM file containing unique valid pairs. */
    @Parameter(names={"-v", "--valid-pairs-bam"}, required = true, description = "Path to BAM file with unique valid pairs produced using the align command.", order = 3)
    private String validPairsBamFile = null;

    /** Path to the genome digest file produced by GOPHER. */
    @Parameter(names={"-d","--digest-file"}, required = true, description = "Path to GOPHER digest file.", order = 4)
    private String digestFile = null;

    /** Aggregate counts for different read pair orienations. */
    @Parameter(names={"-s", "--split-counts"},description = "Split counts for different read pair orientations.", order = 5)
    private boolean split=false;


    public CountCommand() {
    }

    public void execute() throws DiachromaticException {

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
            counter.printFragmentInteractionCountsMapAsCountTable();
            counter.printStatistics();
            logger.trace("...done!");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }
    @Override
    public String toString() {return "diachromatic:count";} //???
}
