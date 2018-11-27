package org.jax.diachromatic.count;

import htsjdk.samtools.BAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.InteractionCountsMap;

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
     * Stores counts about interactions.
     */
    InteractionCountsMap interactionMap;

    public Counter(BAMFileReader validPairsBAM, String outputPathPrefix) {

    }

    /**
     * Print statistics to 'prefix.interaction.stats.txt'.
     */
    public void printStatistics(String statsFile) {

    }


}


