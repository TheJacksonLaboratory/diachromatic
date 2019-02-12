package org.jax.diachromatic.command;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.count.Counter;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.FileNotFoundException;

public class CountCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to BAM file containing unique valid pairs. */
    private String validPairsBamFile = null;

    /** Path to the genome digest file produced by GOPHER. */
    private String digestFile = null;

    private String outputPathPrefix = null;

    private String filenamePrefix;


    public CountCommand(String validPairsBamFile, String digestFile, String outputPathPrefix, String filenamePrefix) {
        this.validPairsBamFile = validPairsBamFile;
        this.digestFile = digestFile;
        this.outputPathPrefix = outputPathPrefix;
        this.filenamePrefix = filenamePrefix;
    }

    public void execute() throws DiachromaticException {
        logger.trace(String.format("About to read digests from %s",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);

        Counter counter = new Counter(validPairsBamFile, digestMap, outputPathPrefix, filenamePrefix);
        try {
            counter.countInteractions();
            counter.printStatistics();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

    }
    @Override
    public String toString() {return "diachromatic:count";} //???
}
