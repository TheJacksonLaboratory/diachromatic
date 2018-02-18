package org.jax.diachromatic.truncation;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.util.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class Truncator {
    private static final Logger logger = LogManager.getLogger();

    private final String outdir;
    private final String fastqFile1;
    private final String fastqFile2;
    private final RestrictionEnzyme renzyme;
    private final String filledEndSequence;

    private final String outputSuffix;

    private String outputFASTQ1, outputFASTQ2;
    /**
     * Read 1 was too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int read1tooShort = 0;
    /**
     * Read 2 was too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int read2tooShort;
    /**
     * Read 1 or read 2 or both were too short after truncation, leading to the removal of the affected read inputSAMfiles.
     */
    private int removedBecauseAtLeastOneReadTooShort;

    private static int LENGTH_THRESHOLD = 20;

    /**
     * @param outdir Output directory
     * @param inputFASTQforward input FASTQ file for forward reads
     * @param inputFASTQreverse input FASTQ file for reverse reads
     * @param re restriction enzyme
     * @param suffix suffix for the output file
     */
    public Truncator(String outdir, String inputFASTQforward, String inputFASTQreverse, RestrictionEnzyme re, String suffix) {
        this.outdir = outdir;
        this.fastqFile1 = inputFASTQforward;
        this.fastqFile2 = inputFASTQreverse;
        this.renzyme = re;
        filledEndSequence = fillEnd(renzyme);
        outputSuffix = suffix;
        makeOutdirectoryIfNeeded();
        createOutputNames();
    }

    /**
     * This form of the constructor is used if we want to use a length threshold different from the default value
     * (see {@link #LENGTH_THRESHOLD}).
     * @param outdir Output directory
     * @param inputFASTQforward input FASTQ file for forward reads
     * @param inputFASTQreverse input FASTQ file for reverse reads
     * @param re restriction enzyme
     * @param suffix suffix for the output file
     * @param threshold length threshold for retaining reads after truncation.
     */
    public Truncator(String outdir, String inputFASTQforward, String inputFASTQreverse, RestrictionEnzyme re, String suffix, int threshold) {
        this(outdir,inputFASTQforward,inputFASTQreverse,re,suffix);
        LENGTH_THRESHOLD=threshold;
    }

    private void makeOutdirectoryIfNeeded() {
        File f = new File(outdir);
        if (f.exists() && f.isFile()) {
            logger.error(String.format("Cannot make output directory called %s because a file of the same name exists", outdir));
        } else if (f.exists()) {
            return; // directory already there
        } else {
            f.mkdir();
        }
    }


    private void createOutputNames() {
        String basename1 = (new File(fastqFile1)).getName();
        int i = basename1.lastIndexOf(".");
        outputFASTQ1 = String.format("%s%s%s.%s%s", outdir, File.separator, basename1.substring(0, i), outputSuffix, basename1.substring(i));
        String basename2 = (new File(fastqFile2)).getName();
        i = basename2.lastIndexOf(".");
        outputFASTQ2 = String.format("%s%s%s.%s%s", outdir, File.separator, basename2.substring(0, i), outputSuffix, basename2.substring(i));
        logger.trace(String.format("F1 %s \n F2 %s", outputFASTQ1, outputFASTQ2));
    }

    /**
     * Parse the two input FASTQ files using a {@link FastqPairParser} object that returns one inputSAMfiles of reads at a time.
     * For each inputSAMfiles of reads, if one or both of the reads was truncated to the extent that the remaining read is too
     * short, then skip the read inputSAMfiles. Write out each read of valid pairs to separate output files.
     */
    public void parseFASTQ() {
        FastQRecord.setLigationSequence(filledEndSequence);
        FastQRecord.setRestrictionSequence(renzyme.getPlainSite());
        FastqPairParser parser = new FastqPairParser(fastqFile1, fastqFile2, filledEndSequence);
        removedBecauseAtLeastOneReadTooShort = 0;
        try {
            BufferedWriter out1 = new BufferedWriter(new FileWriter(outputFASTQ1));
            BufferedWriter out2 = new BufferedWriter(new FileWriter(outputFASTQ2));
            while (parser.hasNextPair()) {
                Pair<FastQRecord, FastQRecord> pair = parser.getNextPair();
                if (pair.first.getLen() < LENGTH_THRESHOLD) {
                    if (pair.second.getLen() < LENGTH_THRESHOLD) read2tooShort++;
                    read1tooShort++;
                    removedBecauseAtLeastOneReadTooShort++;
                    continue;
                } else if (pair.second.getLen() < LENGTH_THRESHOLD) {
                    read2tooShort++;
                    removedBecauseAtLeastOneReadTooShort++;
                } else {
                    pair.first.writeToStream(out1);
                    pair.second.writeToStream(out2);
                }
            }
            out1.close();
            out2.close();
        } catch (IOException e) {
            logger.fatal(String.format("Error encountered while writing truncated FASTQ files: %s", e.getMessage()));
            e.printStackTrace();
        }
        logger.trace(String.format("Number of reads processed: %d and Number of forward reads truncated %d (%.2f%%)",
                parser.getnReadsProcessed(), parser.getReadOneTruncated(),
                100.0 * parser.getReadOneTruncated() / parser.getnReadsProcessed()));
        logger.trace(String.format("removed b/c too short %d", removedBecauseAtLeastOneReadTooShort));

    }

    /**
     * The ligation sequence in capture Hi-C is the result of cutting DNA with a restriction enzyme, filling in the
     * overhands with biotinylated nucleotides, and performing blunet ended ligation. For examples, HindIII has the
     * following cutting sequence: {@code HindIII A^AGCTT}, thus, it cuts between the first and second nucleotides (both A).
     * The ligation seqeunce will be {@code A + AGCT + AGCT + T = AAGCTAGCTT}.
     *
     * @param re restriction enzyme
     * @return
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
