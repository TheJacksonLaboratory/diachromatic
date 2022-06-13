package org.jax.diachromatic.count;


import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;

import org.jax.diachromatic.align.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 * This class is intended for counting read pairs for piars of restriction fragments and for counting reads at
 * interacting fragments. This is currently done in class Align, but should be moved to this class for better
 * clarity.
 *
 * The input for the constructor will be a BAM file containing the valid read pairs as well as a prefix
 * for summarize file (including the path).
 *
 * The summarize will consist of three files:
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
public class BadReadsCounter {
    private static final Logger logger = LoggerFactory.getLogger(Counter.class);
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * HashMap that stores interaction counts. Key: reference of digest pair; Value: SimpleTwistedCount objects.
     */
    private Map<DigestPair,BadReadCount> dp2countsMap;

    /**
     * Stores interaction counts.
     */
    private final DigestMap digestMap;

    /**
     * Paths for summarize files.
     */
    private String outputTsvInteractingFragmentCounts;
    private String outputTsvInteractionCounts;
    private String outputWashUSimpleInteractionCounts;
    private String outputTxtStats;

    /**
     * A reader for the unique valid read pairs.
     */
    final private SamReader reader;

    /**
     * Iterator over reads from {@link #reader}.
     */
    final private Iterator<SAMRecord> it;

    /**
     * Keeps track of the number of processed read pairs
     */
    private int n_pairs_total = 0;

    private int n_unpaired = 0;

    private int n_same_digest = 0;


    /**
     * Total number of interactions
     */
    private int interaction_count = 0;



    /**
     * Total number of active interacting fragments
     */
    private int active_interacting_fragment_count = 0;

    /**
     * Largest number of read pairs for given digest pairs.
     */
    private static int MAX_K = 20000;

    /**
     * Largest number of read pairs for given digest pairs.
     */
    private static int LONG_RANGE_THRESHOLD = 10000;

    /**
     * Array for counting interactions with k read pairs. The index corresponds to k, e.g. array[2]
     * contains the number of interactions with 2 read pairs.
     */
    private int[] kInteractionCounts =  new int[MAX_K+1];

    boolean split = false;


    public BadReadsCounter(SamReader samReader, DigestMap digestMap,  String outputDirAndFilePrefix) {
        this.reader = samReader;
        this.digestMap = digestMap;
        this.it = reader.iterator();
        createOutputNames(outputDirAndFilePrefix);
        this.dp2countsMap = new HashMap<>();
    }

    public void countInteractions() {

        // iterate over unique valid pairs
        n_pairs_total = 0;
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();
            ReadPair readPair = new ReadPair(record1, record2, digestMap);
            n_pairs_total ++;
            //readPair.setRelativeOrientationTag();
            if (! readPair.isPaired()) {
                // at least one read is not uniquely mapped, skip this readpair
                n_unpaired++;
                continue;
            }
            DigestPair dp = readPair.getDigestPair();
            if (dp.isSameDigest()) {
                n_same_digest++;
                continue;
            }
            //record1.get
            //record1.get
            incrementDigestPair(dp, readPair);

            if (interaction_count % 10000000 == 0) {
                logger.trace("Number of Interactions: " + interaction_count);
            }


            n_pairs_total++;
        }
    }

    public void incrementDigestPair(DigestPair dp, ReadPair rp) {

        if (!dp2countsMap.containsKey(dp)) {
            // this is the first read pair for this pair of digests
            interaction_count++;

           // dp2countsMap.put(dp, new SimpleTwistedCount());
        }

    }


    public int getInteractionCount(){
        return interaction_count;
    }

    /**
     * Print statistics to 'prefix.interaction.stats.txt'.
     */
    public void printStatistics() throws FileNotFoundException {

        // create file for summarize
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTxtStats));

        printStream.print("#Count statistics\n");
        printStream.print("==================\n\n");
        printStream.print("total_read_pairs_ processed:" + n_pairs_total + "\n");
        //  Counts of read pair orientations

        // Summary statistics about interactions between active (most typically enriched) and inactive fragments
        printStream.print("total_interaction_count:" + interaction_count + "\n");
           printStream.print("");
        printStream.print("selected_interacting_fragments:" + active_interacting_fragment_count + "\n");
        // Quality metrics for experimental trouble shooting:
       // printStream.print("target_enrichment_coefficient:" + String.format("%.2f%%", 100*this.getTargetEnrichmentCoefficient()) + "\n");
//        double fsi = 100.0*n_singleton_interactions/interaction_count;
//        printStream.print("fraction_singleton_interactions:" + String.format("%.2f%%", fsi) + "\n");
          printStream.print("\n");
        printStream.print("self_ligated_fragment_size_count_array:");
        for(int i=0; i<MAX_K; i++) {
            if (i < MAX_K - 1) {
                printStream.print(kInteractionCounts[i] + ", ");

            } else {
                printStream.print(kInteractionCounts[i]);

            }
        }
        printStream.print("\n");
    }

    private void createOutputNames(String outputPathPrefix) {
        outputTsvInteractingFragmentCounts = String.format("%s.%s", outputPathPrefix, "interacting.fragments.counts.table.tsv");
        outputTsvInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts.table.tsv");
        outputWashUSimpleInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts.washU.simple.tsv");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "count.stats.txt");
    }


    /**
     * Prints digest pairs with associated read pair counts to a tab separated file.
     *
     * @throws FileNotFoundException  if the file output stream cannot be open for the TSV file of interaction counts
     */
    public void printInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // init array for k-interaction counting
        Arrays.fill(kInteractionCounts, 0);

        // create file for summarize
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTsvInteractionCounts));



    }







}

