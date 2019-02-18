package org.jax.diachromatic.count;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.*;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

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
    //private InteractionCountsMap interactionMap;


    private Map<DigestPair,SimpleTwistedCount> dp2countsMap;

    /**
     * Stores interaction counts.
     */
    private DigestMap digestMap;

    /**
     * Paths for output files.
     */
    private String outputTsvInteractingFragmentCounts;
    private String outputTsvInteractionCounts;
    private String outputTsvInteractionCounts2;
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


    public Counter(String validPairsBamFile, DigestMap digestMap, String outputPathPrefix, String outputDirAndFilePrefix) {
        this.reader = SamReaderFactory.makeDefault().open(new File(validPairsBamFile));
        this.digestMap=digestMap;
        this.it = reader.iterator();
        createOutputNames(outputDirAndFilePrefix);
        this.dp2countsMap=new HashMap<>();
    }

    public void countInteractions() throws DiachromaticException, FileNotFoundException {

        //logger.trace("About to determine interaction counts...");
        //interactionMap = new InteractionCountsMap();

        // iterate over unique valid pairs
        n_pairs_total = 0;
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();

            // create read pair
            ReadPair readPair = new ReadPair(record1, record2, digestMap);

            read_count=read_count+2;
            if(readPair.forwardDigestIsActive()){active_read_count++;}
            if(readPair.reverseDigestIsActive()){active_read_count++;}

            // count interaction
            /*interactionMap.incrementFragPair(
                    readPair.forward().getReferenceName(),
                    readPair.getForwardDigestStart(),
                    readPair.getForwardDigestEnd(),
                    readPair.forwardDigestIsActive(),
                    readPair.reverse().getReferenceName(),
                    readPair.getReverseDigestStart(),
                    readPair.getReverseDigestEnd(),
                    readPair.reverseDigestIsActive(),
                    readPair.getRelativeOrientationTag());*/

            DigestPair dp = readPair.getDigestPair();
            if(!dp2countsMap.containsKey(dp)) {
                // this is the first read pair for this pair of digests
                interaction_count++;
                if (readPair.forwardDigestIsActive() && readPair.reverseDigestIsActive()) {
                    active_active_interaction_count++;
                } else if (!readPair.forwardDigestIsActive() && !readPair.reverseDigestIsActive()) {
                    inactive_inactive_interaction_count++;
                } else {
                    active_inactive_interaction_count++;
                }
                dp2countsMap.put(dp,new SimpleTwistedCount());
            }
            if (readPair.isTwisted()) {
                dp2countsMap.get(dp).twisted++;
            } else {
                dp2countsMap.get(dp).simple++;
            }


            if(interaction_count%1000000==0) {
                logger.trace("Number of Interactions: " + interaction_count);
            }

            if(readPair.getRelativeOrientationTag().equals("F1F2")) {n_F1F2++;}
            if(readPair.getRelativeOrientationTag().equals("F2F1")) {n_F2F1++;}
            if(readPair.getRelativeOrientationTag().equals("R1R2")) {n_R1R2++;}
            if(readPair.getRelativeOrientationTag().equals("R2R1")) {n_R2R1++;}
            if(readPair.getRelativeOrientationTag().equals("F1R2")) {n_F1R2++;}
            if(readPair.getRelativeOrientationTag().equals("R2F1")) {n_R2F1++;}
            if(readPair.getRelativeOrientationTag().equals("F2R1")) {n_F2R1++;}
            if(readPair.getRelativeOrientationTag().equals("R1F2")) {n_R1F2++;}

            if(readPair.isTrans()) {
                n_trans_pairs++;
            }

            n_pairs_total++;
        }
        //logger.trace("...done with counting!");
        //logger.trace("About to print the results...");
        //logger.trace("interactionMap.getTotalNumberOfInteractions(): " + interactionMap.getTotalNumberOfInteractions());
        //interactionMap.printInteractionCountsMapAsCountTable(outputTsvInteractionCounts);
        //interactionMap.printFragmentInteractionCountsMapAsCountTable(outputTsvInteractingFragmentCounts);
        //logger.trace("...done!");
    }

    /**
     * Print statistics to 'prefix.interaction.stats.txt'.
     */
    public void printStatistics() throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTxtStats));

        printStream.print("Summary statistics\n");
        printStream.print("==================\n\n");
        printStream.print("\n");
        printStream.print("Interaction count statistics\n");
        printStream.print("----------------------------\n");
        printStream.print("\n");
        printStream.print("Total number of read pairs processed:\t" + n_pairs_total + "\n\n");

        printStream.print("Counts of read pair orientations:\n");
        printStream.print("\tF1F2 - commie:\t" + n_F1F2 + String.format(" (%.2f%%)", 100.0*n_F1F2/n_pairs_total) + "\n");
        printStream.print("\tF2F1 - commie:\t" + n_F2F1 + String.format(" (%.2f%%)", 100.0*n_F2F1/n_pairs_total) + "\n");
        printStream.print("\tR1R2 - commie:\t" + n_R1R2 + String.format(" (%.2f%%)", 100.0*n_R1R2/n_pairs_total) + "\n");
        printStream.print("\tR2R1 - commie:\t" + n_R2R1 + String.format(" (%.2f%%)", 100.0*n_R2R1/n_pairs_total) + "\n");
        printStream.print("\tF1R2 - innie:\t" + n_F1R2 + String.format(" (%.2f%%)", 100.0*n_F1R2/n_pairs_total) + "\n");
        printStream.print("\tF2R1 - innie:\t" + n_F2R1 + String.format(" (%.2f%%)", 100.0*n_F2R1/n_pairs_total) + "\n");
        printStream.print("\tR2F1 - outie:\t" + n_R2F1 + String.format(" (%.2f%%)", 100.0*n_R2F1/n_pairs_total) + "\n");
        printStream.print("\tR1F2 - outie:\t" + n_R1F2 + String.format(" (%.2f%%)", 100.0*n_R1F2/n_pairs_total) + "\n");
        printStream.print("\n");

        printStream.print("\n");
        printStream.print("Summary statistics about interactions between active and inactive fragments:\n");
        printStream.print("\t" + "Total number of interactions: " + interaction_count + "\n");
        printStream.print("\t" + "Number of interactions between active fragments: " + active_active_interaction_count + "\n");
        printStream.print("\t" + "Number of interactions between inactive fragments: " + inactive_inactive_interaction_count + "\n");
        printStream.print("\t" + "Number of interactions between active and inactive fragments: " + active_inactive_interaction_count + "\n");
        printStream.print("");
        printStream.print("\t" + "Total number of interacting fragments: " + interacting_fragment_count + "\n");
        printStream.print("\t" + "Number of active interacting fragments: " + active_interacting_fragment_count + "\n");
        printStream.print("\n");

        printStream.print("Quality metrics for experimental trouble shooting:\n");
        printStream.print("\tTarget Enrichment Coefficient (TEC):\t" + String.format("%.2f%%", 100*this.getTargetEnrichmentCoefficient()) + "\n");
        printStream.print("\tCross-ligation coefficient (CLC):\t" + String.format("%.2f%%", 100.0*n_trans_pairs/n_pairs_total) + "\n");
        double fsi = 100.0*n_singleton_interactions/interaction_count;
        printStream.print("\tFraction of Singleton Interactions (FSI):\t" + String.format("%.2f%%", fsi) + "\n");
    }

    private void createOutputNames(String outputPathPrefix) {
        outputTsvInteractingFragmentCounts = String.format("%s.%s", outputPathPrefix, "interacting.fragments.counts.table.tsv");
        outputTsvInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts.table.tsv");
        outputTsvInteractionCounts2 = String.format("%s.%s", outputPathPrefix, "interaction.counts2.table.tsv");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "count.stats.txt");
    }


    /**
     * Prints digest pairs with associated read pair counts to a tab separated file.
     *
     * @throws FileNotFoundException
     */
    public void printInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTsvInteractionCounts2));

        for (DigestPair dp : this.dp2countsMap.keySet()) {
            SimpleTwistedCount cc = this.dp2countsMap.get(dp);
            printStream.println(dp.toString() + "\t" + cc.simple + ":" + cc.twisted);
            if(cc.simple+cc.twisted==1) {this.n_singleton_interactions++;}
        }
    }

    /**
     * Prints coordinates of interacting digests and associated read counts to a tab separated file.
     *
     * @throws FileNotFoundException
     */
    public void printFragmentInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTsvInteractingFragmentCounts));

        HashMap<Digest,Integer> readCountsAtDigestsMap = new HashMap<>();

        // Iterate over all interactions and add individual digest to a hashMap with key=digestRef and value=read count.
        Integer readCount;
        for (DigestPair dp : this.dp2countsMap.keySet()) {
            SimpleTwistedCount cc = this.dp2countsMap.get(dp);
            if(!readCountsAtDigestsMap.containsKey(dp.forward())) {
                readCountsAtDigestsMap.put(dp.forward(), 1);
                interacting_fragment_count++;
                if(dp.forward().isSelected()) {
                    active_interacting_fragment_count++;
                }
            } else {
                readCount = cc.simple + cc.twisted;
                readCountsAtDigestsMap.put(dp.forward(), readCountsAtDigestsMap.get(dp.forward()) + readCount);
            }
            if(!readCountsAtDigestsMap.containsKey(dp.reverse())) {
                readCountsAtDigestsMap.put(dp.reverse(), 1);
                interacting_fragment_count++;
                if(dp.reverse().isSelected()) {
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
            if(d.isSelected()) {c = 'A';}
            printStream.println(d.getChromosome() + "\t" + d.getDigestStartPosition() + "\t" + d.getDigestEndPosition() + "\t" + c + "\t" + readCountsAtDigestsMap.get(d));
        }

    }


    /**
     * @return Percentage of reads in selective/active digests.
     */
    public double getTargetEnrichmentCoefficient() {
        return 1.0*active_read_count/read_count;
    }



}


