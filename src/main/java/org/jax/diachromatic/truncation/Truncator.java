package org.jax.diachromatic.truncation;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.util.Pair;

import java.io.*;
import java.util.zip.GZIPOutputStream;

public class Truncator {
    private static final Logger logger = LogManager.getLogger();

    private final String fastqFile1;
    private final String fastqFile2;
    private final RestrictionEnzyme renzyme;
    private final String filledEndSequence;

    private String outputFASTQ1, outputFASTQ2, outputSummaryStatistics;
    /**
     * Read 1 was too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int removedBecauseRead1TooShort;
    /**
     * Read 2 was too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int removedBecauseRead2TooShort;
    /**
     * Read 1 or read 2 or both were too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int NumOfPairsRemovedBecauseAtLeastOneReadTooShort;
    /**
     * Count reads that start with the dangling end sequence.
     */
    private int numOfMaybeDanglingRead1;
    private int numOfMaybeDanglingRead2;

    private static final int LENGTH_THRESHOLD = 19; // using 19 the same results as for HiCUP are obtained

    public Truncator(String inputFASTQforward, String inputFASTQreverse, RestrictionEnzyme re, boolean stickyEnds, String outputPathPrefix) {
        this.fastqFile1 = inputFASTQforward;
        this.fastqFile2 = inputFASTQreverse;
        this.renzyme = re;
        if(stickyEnds) {
            filledEndSequence=re.getPlainSite();
        } else {
            filledEndSequence = fillEnd(renzyme);
        }
        outputFASTQ1 = String.format("%s.%s", outputPathPrefix, "truncated_R1.fastq.gz");
        outputFASTQ2 = String.format("%s.%s", outputPathPrefix, "truncated_R2.fastq.gz");
        outputSummaryStatistics = String.format("%s.%s", outputPathPrefix, "truncation.stats.txt");
    }

    /**
     * Parse the two input FASTQ files using a {@link FastqPairParser} object that returns one inputSAMfiles of reads at a time.
     * For each inputSAMfiles of reads, if one or both of the reads was truncated to the extent that the remaining read is too
     * short, then skip the read inputSAMfiles. Write out each read of valid pairs to separate output files.
     */
    public void parseFASTQ() throws DiachromaticException {
        PotentiallyTruncatedFastQRecord.setLigationSequence(filledEndSequence);
        PotentiallyTruncatedFastQRecord.setRestrictionSequence(renzyme.getPlainSite());
        PotentiallyTruncatedFastQRecord.setDanglingSequence(renzyme.getDanglingEndSequence());
        FastqPairParser parser = new FastqPairParser(fastqFile1, fastqFile2, filledEndSequence);
        NumOfPairsRemovedBecauseAtLeastOneReadTooShort = 0;
        removedBecauseRead1TooShort = 0;
        removedBecauseRead2TooShort = 0;
        numOfMaybeDanglingRead1 = 0;
        numOfMaybeDanglingRead2 = 0;
        logger.trace("filledEndSequence:"  + filledEndSequence + "\trenzyme.getSite(): " + renzyme.getSite() + "\tenzyme.getPlainSite(): " + renzyme.getPlainSite() + "\trenzyme.getDanglingEndSequence(): " + renzyme.getDanglingEndSequence() + "\n");
        try {

            BufferedWriter out1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFASTQ1))));
            BufferedWriter out2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFASTQ2))));

            while (parser.hasNextPair()) {
                Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> pair = parser.getNextPair();

                if (pair.first.getLen() < LENGTH_THRESHOLD) {
                    removedBecauseRead1TooShort++;
                }
                if (pair.second.getLen() < LENGTH_THRESHOLD) {
                    removedBecauseRead2TooShort++;
                }
                if((LENGTH_THRESHOLD) < pair.first.getLen() && (LENGTH_THRESHOLD < pair.second.getLen())) {
                    pair.first.writeToStream(out1);
                    pair.second.writeToStream(out2);
                }
                else {
                    NumOfPairsRemovedBecauseAtLeastOneReadTooShort++;
                }
                if(pair.first.isMaybeDangling()) {
                    numOfMaybeDanglingRead1++;
                }
                if(pair.second.isMaybeDangling()) {
                    numOfMaybeDanglingRead2++;
                }
            }
            out1.close();
            out2.close();
        } catch (IOException e) {
            logger.fatal(String.format("Error encountered while writing truncated FASTQ files: %s", e.getMessage()));
            e.printStackTrace();
        }
        logger.trace(String.format("Number of pairs processed: %d",
                parser.getnReadsProcessed()));
        logger.trace(String.format("Number of truncated forward reads: %d (%.2f%%)",
                parser.getReadOneTruncated(),
                100.0 * parser.getReadOneTruncated() / parser.getnReadsProcessed()));
        logger.trace(String.format("Number of truncated reverse reads: %d (%.2f%%)",
                parser.getReadTwoTruncated(),
                100.0 * parser.getReadOneTruncated() / parser.getnReadsProcessed()));
        logger.trace(String.format("Number of maybe dangling forward reads: %d (%.2f%%)", numOfMaybeDanglingRead1,100.0 * numOfMaybeDanglingRead1/parser.getnReadsProcessed()));
        logger.trace(String.format("Number of maybe dangling reverse reads: %d (%.2f%%)", numOfMaybeDanglingRead2,100.0 * numOfMaybeDanglingRead2/parser.getnReadsProcessed()));
        logger.trace(String.format("Number of too short removed forward reads (<%d): %d", LENGTH_THRESHOLD, removedBecauseRead1TooShort));
        logger.trace(String.format("Number of too short removed reverse reads (<%d): %d", LENGTH_THRESHOLD, removedBecauseRead2TooShort));
        logger.trace(String.format("Number of removed pairs (at least one read too short): %d (%.2f%%)",
                NumOfPairsRemovedBecauseAtLeastOneReadTooShort,
                100.0 * NumOfPairsRemovedBecauseAtLeastOneReadTooShort / parser.getnReadsProcessed()));

        PrintStream printSummaryStatistics = null;
        try {
            printSummaryStatistics = new PrintStream(new FileOutputStream(outputSummaryStatistics));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        printSummaryStatistics.print("Summary statistics\n");
        printSummaryStatistics.print("==================\n\n");
        printSummaryStatistics.print("\n");
        printSummaryStatistics.print("Truncation statistics\n");
        printSummaryStatistics.print("---------------------\n");
        printSummaryStatistics.print("\n");
        printSummaryStatistics.print("Total number of read pairs processed:\t" + parser.getnReadsProcessed() + "\n\n");
        printSummaryStatistics.print(String.format("Number of truncated forward reads: %d (%.2f%%)\n",
                parser.getReadOneTruncated(),
                100.0 * parser.getReadOneTruncated() / parser.getnReadsProcessed()));
        printSummaryStatistics.print(String.format("Number of truncated reverse reads: %d (%.2f%%)\n\n",
                parser.getReadTwoTruncated(),
                100.0 * parser.getReadOneTruncated() / parser.getnReadsProcessed()));
        printSummaryStatistics.print(String.format("Number of maybe dangling forward reads: %d (%.2f%%)\n", numOfMaybeDanglingRead1,100.0 * numOfMaybeDanglingRead1/parser.getnReadsProcessed()));
        printSummaryStatistics.print(String.format("Number of maybe dangling reverse reads: %d (%.2f%%)\n\n", numOfMaybeDanglingRead2,100.0 * numOfMaybeDanglingRead2/parser.getnReadsProcessed()));
        printSummaryStatistics.print(String.format("Number of too short removed forward reads (<%d): %d\n", LENGTH_THRESHOLD, removedBecauseRead1TooShort));
        printSummaryStatistics.print(String.format("Number of too short removed reverse reads (<%d): %d\n", LENGTH_THRESHOLD, removedBecauseRead2TooShort));
        printSummaryStatistics.print(String.format("Number of removed pairs (at least one read too short): %d (%.2f%%)",
                NumOfPairsRemovedBecauseAtLeastOneReadTooShort,
                100.0 * NumOfPairsRemovedBecauseAtLeastOneReadTooShort / parser.getnReadsProcessed()));
    }


    /**
     * The ligation sequence in capture Hi-C is the result of cutting DNA with a restriction enzyme, filling in the
     * overhands with biotinylated nucleotides, and performing blunet ended ligation. For examples, HindIII has the
     * following cutting sequence: {@code HindIII A^AGCTT}, thus, it cuts between the first and second nucleotides (both A).
     * The ligation seqeunce will be {@code A + AGCT + AGCT + T = AAGCTAGCTT}.
     *
     * @param re restriction enzyme
     * @return a string representing the filled-in ligation sequence
     */
    static String fillEnd(RestrictionEnzyme re) {
        String plainsite = re.getPlainSite();
        int offset = re.getOffset();
        int len = plainsite.length();

        if (offset == 0) {
            // this means the enzyme cuts right before the recognition site, e.g., DpnII (5'-^GATC-3')
            return String.format("%s%s", plainsite, plainsite);

        } else if (offset == len) {
            // this means the enzyme cuts right after the recognition site, e.g., NlaIII (5'-CATG^-3')
            return String.format("%s%s", plainsite, plainsite);
        } else {
            int number_of_bases_before_cut = offset;
            int number_of_bases_after_cut = len - offset;
            int flank_size;
            if (number_of_bases_before_cut > len / 2) {
                flank_size = number_of_bases_after_cut;
            } else {
                flank_size = number_of_bases_before_cut;
            }
            // logger.trace(String.format("cut %s n before cut %d and after cut %s  flank size=%d",re.getName(),number_of_bases_before_cut,number_of_bases_after_cut,flank_size ));
            String fivePrimeFlank = plainsite.substring(0, flank_size);
            String fillIn = plainsite.substring(flank_size, len - flank_size); // sequence created by filling in the overhang
            String threePrimeFlank = plainsite.substring(len - flank_size);

            return String.format("%s%s%s%s", fivePrimeFlank, fillIn, fillIn, threePrimeFlank);
        }
    }
}
