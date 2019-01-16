package org.jax.diachromatic.count;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.align.InteractionCountsMap;
import org.jax.diachromatic.align.ReadPair;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;
import java.util.Iterator;

/**
 * This class is intended for counting read pairs for piars of restriction fragments and for counting reads at
 * interacting fragments. This is currently done in class Align, but should be moved to this class for better
 * clarity.
 *
 * The input for the constructor will be a BAM file containing the valid read pairs as well as a prefix
 * for output file (including the path).
 *
 * The output will consist of three files:
 *
 * <li>prefix.interacting.fragments.counts.table.tsv</li>
 * <li>prefix.interaction.counts.table.tsv</li>
 * <li>prefix.interaction.stats.txt</li>
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.0 (2018-05-03)
 *
 */
public class Counter {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * Stores interaction counts.
     */
    InteractionCountsMap interactionMap;

    /**
     * Stores interaction counts.
     */
    DigestMap digestMap;

    /**
     * Paths for output files.
     */
    String outputTsvInteractingFragmentCounts, outputTsvInteractionCounts, outputTxtStats;

    /**
     * A reader for the unique valid read pairs.
     */
    final private SamReader reader;

    /**
     * Iterator over reads from {@link #reader}.
     */
    final private Iterator<SAMRecord> it;


    public Counter(String validPairsBamFile, DigestMap digestMap, String outputPathPrefix, String filenamePrefix) {
        this.reader = SamReaderFactory.makeDefault().open(new File(validPairsBamFile));
        this.digestMap=digestMap;
        this.it = reader.iterator();

    }

    public void countInteractions() throws DiachromaticException {

        interactionMap = new InteractionCountsMap(1);

        // iterate over unique valid pairs
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();
            logger.trace(record1.getReadName());
            logger.trace(record1.getAttribute("RO"));

            // create read pair
            ReadPair readPair = new ReadPair(record1, record2, digestMap);
            //readPair.

        }

    }

    /**
     * Print statistics to 'prefix.interaction.stats.txt'.
     */
    public void printStatistics() {

    }

    private void createOutputNames(String outputPathPrefix) {
        outputTsvInteractingFragmentCounts = String.format("%s.%s", outputPathPrefix, "interacting.fragments.counts.table.tsv"); // will be moved to class counts
        outputTsvInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts.table.tsv"); // will be moved to class counts
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "count.stats.txt");
    }



}


