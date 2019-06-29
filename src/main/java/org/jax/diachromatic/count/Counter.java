package org.jax.diachromatic.count;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.*;

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
public class Counter {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * HashMap that stores interaction counts. Key: reference of digest pair; Value: SimpleTwistedCount objects.
     */
    private Map<DigestPair,SimpleTwistedCount> dp2countsMap;

    /**
     * Stores interaction counts.
     */
    private DigestMap digestMap;

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
    private int n_pairs_total;

    /**
     * Keeps track of the number of processed trans read pairs
     */
    private int  n_trans_pairs;

    /**
     * Count variables for different orientations of read pairs
     */
    private int n_F1F2 = 0;
    private int n_F2F1 = 0;
    private int n_R1R2 = 0;
    private int n_R2R1 = 0;
    private int n_F1R2 = 0;
    private int n_R1F2 = 0;
    private int n_R2F1 = 0;
    private int n_F2R1 = 0;

    /**
     * Total number of reads
     */
    private int read_count = 0;

    /**
     * Total number of reads within active fragments
     */
    private int active_read_count = 0;

    /**
     * Total number of interactions
     */
    private int interaction_count = 0;


    /**
     * Total current number of singleton interactions
     * Note: Initialized only after execution of the function 'printInteractionCountsMapAsCountTable'.
     */
    private int n_singleton_interactions = 0;
    private int n_singleton_interactions_trans = 0;
    private int n_singleton_interactions_short_range = 0;
    private int n_singleton_interactions_long_range = 0;

    private int n_gt1_interaction_count = 0;
    private int n_gt1_interaction_count_trans = 0;
    private int n_gt1_interaction_count_short_range = 0;
    private int n_gt1_interaction_count_long_range = 0;

    /**
     * Number of interactions between two active fragments
     */
    private int active_active_interaction_count = 0;

    /**
     * Number of interactions between two inactive fragments
     */
    private int inactive_inactive_interaction_count = 0;

    /**
     * Number of interactions between active and inactive fragments (both directions)
     */
    private int active_inactive_interaction_count = 0;

    /**
     * Total number of interacting fragments
     */
    private int interacting_fragment_count = 0;

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


    public Counter(SamReader samReader, DigestMap digestMap, String outputDirAndFilePrefix, boolean split) {
        this.reader = samReader;
        this.digestMap = digestMap;
        this.it = reader.iterator();
        createOutputNames(outputDirAndFilePrefix);
        this.dp2countsMap = new HashMap<>();
        this.split=split;
    }

    public void countInteractions() {

        // iterate over unique valid pairs
        n_pairs_total = 0;
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();

            // create read pair
            ReadPair readPair = new ReadPair(record1, record2, digestMap);
            //readPair.setRandomRelativeOrientationTag();
            readPair.setRelativeOrientationTag();

            read_count = read_count + 2;
            if (readPair.forwardDigestIsActive()) {
                active_read_count++;
            }
            if (readPair.reverseDigestIsActive()) {
                active_read_count++;
            }

            DigestPair dp = readPair.getDigestPair();
            incrementDigestPair(dp, readPair);

            if (interaction_count % 10000000 == 0) {
                logger.trace("Number of Interactions: " + interaction_count);
            }

            if (readPair.getRelativeOrientationTag().equals("F1F2")) {
                n_F1F2++;
            }
            if (readPair.getRelativeOrientationTag().equals("F2F1")) {
                n_F2F1++;
            }
            if (readPair.getRelativeOrientationTag().equals("R1R2")) {
                n_R1R2++;
            }
            if (readPair.getRelativeOrientationTag().equals("R2R1")) {
                n_R2R1++;
            }
            if (readPair.getRelativeOrientationTag().equals("F1R2")) {
                n_F1R2++;
            }
            if (readPair.getRelativeOrientationTag().equals("R2F1")) {
                n_R2F1++;
            }
            if (readPair.getRelativeOrientationTag().equals("F2R1")) {
                n_F2R1++;
            }
            if (readPair.getRelativeOrientationTag().equals("R1F2")) {
                n_R1F2++;
            }

            if (readPair.isTrans()) {
                n_trans_pairs++;
            }
            n_pairs_total++;
        }
    }

    public void incrementDigestPair(DigestPair dp, ReadPair rp) {

        if (!dp2countsMap.containsKey(dp)) {
            // this is the first read pair for this pair of digests
            interaction_count++;
            if (dp.forward().isSelected() && dp.reverse().isSelected()) {
                active_active_interaction_count++;
            } else if (!dp.forward().isSelected() && !dp.reverse().isSelected()) {
                inactive_inactive_interaction_count++;
            } else {
                active_inactive_interaction_count++;
            }
            dp2countsMap.put(dp, new SimpleTwistedCount());
        }
        if (rp.isTwisted()) {
            dp2countsMap.get(dp).twisted++;
        } else {
            dp2countsMap.get(dp).simple++;
        }
    }

    public SimpleTwistedCount getSimpleTwistedCountForDigestPair(DigestPair dp) {
        return dp2countsMap.get(dp);
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
        printStream.print("\tF1F2_commie:" + n_F1F2 + String.format(" (%.2f%%)", 100.0*n_F1F2/n_pairs_total) + "\n");
        printStream.print("\tF2F1_commie:" + n_F2F1 + String.format(" (%.2f%%)", 100.0*n_F2F1/n_pairs_total) + "\n");
        printStream.print("\tR1R2_commie:" + n_R1R2 + String.format(" (%.2f%%)", 100.0*n_R1R2/n_pairs_total) + "\n");
        printStream.print("\tR2R1_commie:" + n_R2R1 + String.format(" (%.2f%%)", 100.0*n_R2R1/n_pairs_total) + "\n");
        printStream.print("\tF1R2_innie:" + n_F1R2 + String.format(" (%.2f%%)", 100.0*n_F1R2/n_pairs_total) + "\n");
        printStream.print("\tF2R1_innie:" + n_F2R1 + String.format(" (%.2f%%)", 100.0*n_F2R1/n_pairs_total) + "\n");
        printStream.print("\tR2F1_outie:" + n_R2F1 + String.format(" (%.2f%%)", 100.0*n_R2F1/n_pairs_total) + "\n");
        printStream.print("\tR1F2_outie:" + n_R1F2 + String.format(" (%.2f%%)", 100.0*n_R1F2/n_pairs_total) + "\n");
        // Summary statistics about interactions between active (most typically enriched) and inactive fragments
        printStream.print("total_interaction_count:" + interaction_count + "\n");
        printStream.print("interactions_between_selected_fragments:" + active_active_interaction_count + "\n");
        printStream.print("interactions_between_unselected_fragments:" + inactive_inactive_interaction_count + "\n");
        printStream.print("interactions_between_selected_and_unselected_fragments:" + active_inactive_interaction_count + "\n");
        printStream.print("");
        printStream.print("total_interacting_fragments:" + interacting_fragment_count + "\n");
        printStream.print("selected_interacting_fragments:" + active_interacting_fragment_count + "\n");
        // Quality metrics for experimental trouble shooting:
        printStream.print("target_enrichment_coefficient:" + String.format("%.2f%%", 100*this.getTargetEnrichmentCoefficient()) + "\n");
        printStream.print("cross_ligation_coefficient:" + String.format("%.2f%%", 100.0*n_trans_pairs/n_pairs_total) + "\n");
        double fsi = 100.0*n_singleton_interactions/interaction_count;
        printStream.print("fraction_singleton_interactions:" + String.format("%.2f%%", fsi) + "\n");
        printStream.print("n_singleton_interactions:" + n_singleton_interactions + "\n");
        printStream.print("n_singleton_interactions_trans:" + n_singleton_interactions_trans + "\n");
        printStream.print("n_singleton_interactions_short_range:" + n_singleton_interactions_short_range + "\n");
        printStream.print("n_singleton_interactions_long_range:" + n_singleton_interactions_long_range + "\n");
        printStream.print("n_gt1_interaction_count:" + n_gt1_interaction_count + "\n");
        printStream.print("n_gt1_interaction_count_trans:" + n_gt1_interaction_count_trans + "\n");
        printStream.print("n_gt1_interaction_count_short_range:" + n_gt1_interaction_count_short_range + "\n");
        printStream.print("n_gt1_interaction_count_long_range:" + n_gt1_interaction_count_long_range + "\n");
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

        for (DigestPair dp : this.dp2countsMap.keySet()) {
            SimpleTwistedCount cc = this.dp2countsMap.get(dp);
            kInteractionCounts[cc.simple + cc.twisted]++;
            //int cnt = cc.simple + cc.twisted;
            //printStream.println(dp.toString() + "\t" + cnt);
            if(this.split) {
                printStream.println(dp.toString() + "\t" + cc.simple + ":" + cc.twisted);
            } else {
                int c = cc.simple + cc.twisted;
                printStream.println(dp.toString() + "\t" + c);
            }
            if (cc.simple + cc.twisted == 1) {
                this.n_singleton_interactions++;
                if (!dp.forward().getChromosome().equals(dp.reverse().getChromosome())) {
                    n_singleton_interactions_trans++;
                } else {
                    int forward_digest_center = dp.forward().getDigestStartPosition() + ((dp.forward().getDigestEndPosition() - dp.forward().getDigestStartPosition()) / 2);
                    int reverse_digest_center = dp.reverse().getDigestStartPosition() + ((dp.reverse().getDigestEndPosition() - dp.reverse().getDigestStartPosition()) / 2);
                    if(Math.abs(reverse_digest_center - forward_digest_center)<LONG_RANGE_THRESHOLD) {
                        n_singleton_interactions_short_range++;
                    } else {
                        n_singleton_interactions_long_range++;
                    }
                }
            } else {
                n_gt1_interaction_count++;
                if (!dp.forward().getChromosome().equals(dp.reverse().getChromosome())) {
                    n_gt1_interaction_count_trans++;
                } else {
                    int forward_digest_center = dp.forward().getDigestStartPosition() + ((dp.forward().getDigestEndPosition() - dp.forward().getDigestStartPosition()) / 2);
                    int reverse_digest_center = dp.reverse().getDigestStartPosition() + ((dp.reverse().getDigestEndPosition() - dp.reverse().getDigestStartPosition()) / 2);
                    if(Math.abs(reverse_digest_center - forward_digest_center)<LONG_RANGE_THRESHOLD) {
                        n_gt1_interaction_count_short_range++;
                    } else {
                        n_gt1_interaction_count_long_range++;
                    }
                }
            }
        }
    }

    /**
     * Prints digest pairs with associated read pair counts to simple text format established by washU.
     *
     * @throws FileNotFoundException  if the file output stream cannot be opened.
     */
    public void printInteractionCountsMapInWashUSimpleTextFormat() throws FileNotFoundException {

        // create file
        PrintStream printStream = new PrintStream(new FileOutputStream(outputWashUSimpleInteractionCounts));

        for (DigestPair dp : this.dp2countsMap.keySet()) {
            if(dp.forward().getChromosome().equals(dp.reverse().getChromosome())){
                int forward_digest_center = dp.forward().getDigestStartPosition() + ((dp.forward().getDigestEndPosition() - dp.forward().getDigestStartPosition()) / 2);
                int reverse_digest_center = dp.reverse().getDigestStartPosition() + ((dp.reverse().getDigestEndPosition() - dp.reverse().getDigestStartPosition()) / 2);
                if(LONG_RANGE_THRESHOLD<=Math.abs(reverse_digest_center - forward_digest_center)) {
                    SimpleTwistedCount cc = this.dp2countsMap.get(dp);
                    int c = cc.simple + cc.twisted;
                    String coordinatesF = String.format("%s:%s-%s", dp.forward().getChromosome(),dp.forward().getDigestStartPosition(),dp.forward().getDigestEndPosition());
                    String coordinatesR = String.format("%s:%s-%s", dp.reverse().getChromosome(),dp.reverse().getDigestStartPosition(),dp.reverse().getDigestEndPosition());
                    String coordinates;
                    // fragment with the smaller starting position comes first
                    if(dp.forward().getDigestStartPosition()<dp.reverse().getDigestStartPosition()) {
                        coordinates = String.format("%s\t%s", coordinatesF, coordinatesR);
                    } else {
                        coordinates = String.format("%s\t%s", coordinatesR, coordinatesF);
                    }
                    printStream.println(coordinates + "\t" + c);
                }
            }
        }
    }

    /**
     * Prints coordinates of interacting digests and associated read counts to a tab separated file.
     *
     * @throws FileNotFoundException if the file output stream cannot be open for the TSV file of fragment interaction counts
     */
    public void printFragmentInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // create file for summarize
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTsvInteractingFragmentCounts));

        HashMap<Digest, Integer> readCountsAtDigestsMap = new HashMap<>();

        // Iterate over all interactions and add individual digest to a hashMap with key=digestRef and value=read count.
        Integer readCount;
        for (DigestPair dp : this.dp2countsMap.keySet()) {
            SimpleTwistedCount cc = this.dp2countsMap.get(dp);
            if (!readCountsAtDigestsMap.containsKey(dp.forward())) {
                readCountsAtDigestsMap.put(dp.forward(), 1);
                interacting_fragment_count++;
                if (dp.forward().isSelected()) {
                    active_interacting_fragment_count++;
                }
            } else {
                readCount = cc.simple + cc.twisted;
                readCountsAtDigestsMap.put(dp.forward(), readCountsAtDigestsMap.get(dp.forward()) + readCount);
            }
            if (!readCountsAtDigestsMap.containsKey(dp.reverse())) {
                readCountsAtDigestsMap.put(dp.reverse(), 1);
                interacting_fragment_count++;
                if (dp.reverse().isSelected()) {
                    active_interacting_fragment_count++;
                }
            } else {
                readCount = cc.simple + cc.twisted;
                readCountsAtDigestsMap.put(dp.reverse(), readCountsAtDigestsMap.get(dp.reverse()) + readCount);
            }
        }

        // Print unique interacting digests and associated read counts
        for (Digest d : readCountsAtDigestsMap.keySet()) {
            char c = 'I';
            if (d.isSelected()) {
                c = 'A';
            }
            printStream.println(d.getChromosome() + "\t" + d.getDigestStartPosition() + "\t" + d.getDigestEndPosition() + "\t" + c + "\t" + readCountsAtDigestsMap.get(d));
        }

    }


    /**
     * @return Percentage of reads in selective/active digests.
     */
    private double getTargetEnrichmentCoefficient() {
        return 1.0*active_read_count/read_count;
    }



}

