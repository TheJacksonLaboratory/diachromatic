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
@Parameters(commandDescription = "count TODO-more text")
public class CountCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to BAM file containing unique valid pairs. */
    @Parameter(names={"-v", "--valid-pairs-bam"},required = true,description = "path to BAM file with unique valid pairs produced using the 'align' subcommand")
    private String validPairsBamFile = null;

    /** Path to the genome digest file produced by GOPHER. */
    @Parameter(names={"-d","--digest-file"}, required = true,description = "path to GOPHER digest file")
    private String digestFile = null;


    public CountCommand() {
    }

    public void execute() throws DiachromaticException {

        makeOutdirectoryIfNeeded();

        logger.trace(String.format("About to read digests from %s",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);

        String outputDirAndFilePrefix=String.format("%s%s%s", outputDir, File.separator,filenamePrefix);

        SamReader reader = SamReaderFactory.makeDefault().open(new File(validPairsBamFile));

        Counter counter = new Counter(reader, digestMap, outputDirAndFilePrefix);
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
