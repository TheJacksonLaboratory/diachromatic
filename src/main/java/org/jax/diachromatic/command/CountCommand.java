package org.jax.diachromatic.command;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.Digest;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.Counter;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.IOException;
import java.util.List;
import java.util.Map;

public class CountCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to BAM file containing unique valid pairs. */
    private String validPairsBamFile = null;

    /** Path to the genome digest file produced by GOPHER. */
    private String digestFile = null;

    /** Path to BED file containing the coordinates of active digests. */
    private String activeDigestsFile = null;

    private String outputPathPrefix = null;

    private String filenamePrefix;


    public CountCommand(String validPairsBamFile, String digestFile, String activeDigestsFile, String outputPathPrefix, String filenamePrefix) {
        this.validPairsBamFile = validPairsBamFile;
        this.digestFile = digestFile;
        this.activeDigestsFile = activeDigestsFile;
        this.outputPathPrefix = outputPathPrefix;
        this.filenamePrefix = filenamePrefix;
    }

    public void execute() throws DiachromaticException {
        logger.trace(String.format("About to read digests from %s",digestFile));
        DigestMap digestMap = new DigestMap(digestFile, activeDigestsFile);

        Counter counter = new Counter(validPairsBamFile, digestMap, outputPathPrefix, filenamePrefix);
        counter.countInteractions();
        counter.printStatistics();

    }
}
