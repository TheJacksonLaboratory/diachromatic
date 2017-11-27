package org.jax.diachromatic.truncation;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.util.Pair;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public class Truncator  {
    private static final Logger logger = LogManager.getLogger();

    private final String outdir;
    private final String fastqFile1;
    private final String fastqFile2;
    private final RestrictionEnzyme renzyme;
    private final String filledEndSequence;

    private int removedBecauseAtLeastOneReadTooShort;

    private static final int LENGTH_THRESHOLD=20;


    public Truncator(String outdir, String file1, String file2,RestrictionEnzyme re){
        this.outdir=outdir;
        this.fastqFile1=file1;
        this.fastqFile2=file2;
        this.renzyme=re;
        filledEndSequence=fillEnd(renzyme);
    }





    public void parseFASTQ() {
        FastQRecord.setLigationSequence(filledEndSequence);
        FastQRecord.setRestrictionSequence(renzyme.getPlainSite());
        FastqPairParser parser = new FastqPairParser(fastqFile1,fastqFile2,filledEndSequence);
        removedBecauseAtLeastOneReadTooShort=0;
        try {
            BufferedWriter out1=new BufferedWriter(new FileWriter("seq1.fastq"));
            BufferedWriter out2=new BufferedWriter(new FileWriter("seq2.fastq"));
            while (parser.hasNextPair()) {
                Pair<FastQRecord, FastQRecord> pair = parser.getNextPair();
                if (pair.first.getLen()<LENGTH_THRESHOLD || pair.second.getLen()<LENGTH_THRESHOLD) {
                    removedBecauseAtLeastOneReadTooShort++;
                    continue;
                }
                pair.first.writeToStream(out1);
                pair.second.writeToStream(out2);
            }
            out1.close();
            out2.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        logger.trace(String.format("No. reads processed: %d and No. of forward reads truncated %d (%.2f%%)",
                parser.getnReadsProcessed(), parser.getReadOneTruncated(),
                100.0*parser.getReadOneTruncated()/parser.getnReadsProcessed()));
        logger.trace(String.format("removed b/c too short %d",removedBecauseAtLeastOneReadTooShort ));

    }


    /**
     * Number of bases in which the ligation junction and reference genome correspond before first mismatch
     * @param ligationSequence
     * @param re
     * @return

    public static int lengthSame(String ligationSequence, RestrictionEnzyme re) {

    }  */


    /**
     * The ligation sequence in capture Hi-C is the result of cutting DNA with a restriction enzyme, filling in the
     * overhands with biotinylated nucleotides, and performing blunet ended ligation. For examples, HindIII has the
     * following cutting sequence: {@code HindIII A^AGCTT}, thus, it cuts between the first and second nucleotides (both A).
     * The ligation seqeunce will be {@code A + AGCT + AGCT + T = AAGCTAGCTT}.
     * @param re
     * @return
     */
    public static String fillEnd(RestrictionEnzyme re) {
        String plainsite = re.getPlainSite();
        int offset = re.getOffset();
        int len= plainsite.length();

        if (offset==0) {
            // this means the enzyme cuts right before the recognition site, e.g., DpnII (5'-^GATC-3')
            return String.format("%s%s",plainsite,plainsite);
        } else if (offset==len) {
            // this means the enzyme cuts right after the recognition site, e.g., NlaIII (5'-CATG^-3')
            return String.format("%s%s",plainsite,plainsite);
        } else {
            int number_of_bases_before_cut=offset;
            int number_of_bases_after_cut=len-offset;
            int flank_size;
            if (number_of_bases_before_cut > len/2) {
                flank_size= number_of_bases_after_cut;
            } else {
                flank_size=number_of_bases_before_cut;
            }
            // logger.trace(String.format("cut %s n before cut %d and after cut %s  flank size=%d",re.getName(),number_of_bases_before_cut,number_of_bases_after_cut,flank_size ));
            String fivePrimeFlank=plainsite.substring(0,flank_size);
            String fillIn = plainsite.substring(flank_size, len - flank_size); // sequence created by filling in the overhang
            String threePrimeFlank=plainsite.substring(len-flank_size);

            return String.format("%s%s%s%s", fivePrimeFlank, fillIn, fillIn, threePrimeFlank);
        }

    }

}
