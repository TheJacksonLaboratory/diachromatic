package org.jax.diachromatic.count;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.align.ReadPair;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
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
    private InteractionCountsMap interactionMap;

    /**
     * Stores interaction counts.
     */
    private DigestMap digestMap;

    /**
     * Paths for output files.
     */
    private String outputTsvInteractingFragmentCounts;
    private String outputTsvInteractionCounts;
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


    public Counter(String validPairsBamFile, DigestMap digestMap, String outputPathPrefix, String filenamePrefix) {
        this.reader = SamReaderFactory.makeDefault().open(new File(validPairsBamFile));
        this.digestMap=digestMap;
        this.it = reader.iterator();
        createOutputNames(outputPathPrefix);

    }

    public void countInteractions() throws DiachromaticException, FileNotFoundException {

        logger.trace("About to determine interaction counts...");
        interactionMap = new InteractionCountsMap(1);

        // iterate over unique valid pairs
        n_pairs_total = 0;
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();

            // create read pair
            ReadPair readPair = new ReadPair(record1, record2, digestMap);

            // count interaction
            interactionMap.incrementFragPair(0,
                    readPair.forward().getReferenceName(),
                    readPair.getForwardDigestStart(),
                    readPair.getForwardDigestEnd(),
                    readPair.forwardDigestIsActive(),
                    readPair.reverse().getReferenceName(),
                    readPair.getReverseDigestStart(),
                    readPair.getReverseDigestEnd(),
                    readPair.reverseDigestIsActive(),
                    readPair.getRelativeOrientationTag());

            if(interactionMap.getTotalNumberOfInteractionsForCondition(0)%1000000==0) {
                logger.trace("Number of Interactions: " + interactionMap.getTotalNumberOfInteractionsForCondition(0));
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
        logger.trace("...done with counting!");
        logger.trace("About to print the results...");
        logger.trace("interactionMap.getTotalNumberOfInteractionsForCondition(0): " + interactionMap.getTotalNumberOfInteractionsForCondition(0));
        interactionMap.printInteractionCountsMapAsCountTable(outputTsvInteractionCounts);
        //interactionMap.printFragmentInteractionCountsMapAsCountTable(outputTsvInteractingFragmentCounts);
        logger.trace("...done!");
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
        printStream.print("\t" + "Total number of interactions: " + interactionMap.getTotalNumberOfInteractionsForCondition(0) + "\n");
        printStream.print("\t" + "Number of interactions between active fragments: " + interactionMap.getNumberOfInteractionsBetweenActiveFragmentsForCondition(0) + "\n");
        printStream.print("\t" + "Number of interactions between inactive fragments: " + interactionMap.getNumberOfInteractionsBetweenInactiveFragmentsForCondition(0) + "\n");
        printStream.print("\t" + "Number of interactions between active and inactive fragments: " + interactionMap.getNumberOfInteractionsBetweenActiveAndInactiveFragmentsForCondition(0) + "\n");
        printStream.print("");
        printStream.print("\t" + "Total number of interacting fragments: " + interactionMap.getTotalNumberOfInteractingFragmentsForCondition(0) + "\n");
        printStream.print("\t" + "Number of active interacting fragments: " + interactionMap.getTotalNumberOfActiveInteractingFragmentsForCondition(0) + "\n");
        printStream.print("\n");

        printStream.print("Quality metrics for experimental trouble shooting:\n");
        printStream.print("\tTarget Enrichment Coefficient (TEC):\t" + String.format("%.2f%%", 100*interactionMap.getTargetEnrichmentCoefficientForCondition(0)) + "\n");
        printStream.print("\tCross-ligation coefficient (CLC):\t" + String.format("%.2f%%", 100.0*n_trans_pairs/n_pairs_total) + "\n");
        double fsi = 100.0*interactionMap.getTotalNumberOfSingletonInteractionsForCondition(0)/interactionMap.getTotalNumberOfInteractionsForCondition(0);
        printStream.print("\tFraction of Singleton Interactions (FSI):\t" + String.format("%.2f%%", fsi) + "\n");
    }

    private void createOutputNames(String outputPathPrefix) {
        outputTsvInteractingFragmentCounts = String.format("%s.%s", outputPathPrefix, "interacting.fragments.counts2.table.tsv");
        outputTsvInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts2.table.tsv");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "count.stats.txt");
    }



}


