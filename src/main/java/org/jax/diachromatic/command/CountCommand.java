package org.jax.diachromatic.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.Counter;
import org.jax.diachromatic.exception.DiachromaticException;

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
        logger.trace(String.format("About to read digests from %s",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);

        Counter counter = new Counter(validPairsBamFile, digestMap, outputPath, filenamePrefix);
        try {
            counter.countInteractions();
            counter.printStatistics();
            counter.printDigestPairsNEW();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }
    @Override
    public String toString() {return "diachromatic:count";} //???
}
