package org.jax.diachromatic.truncation;

import htsjdk.samtools.fastq.FastqRecord;

import java.io.BufferedWriter;
import java.io.IOException;

/**
 * This class is initialized by a HTSJDK FastqRecord object. These obejcts are immutable, and difficult to If the sequence contains a ligation sequence, it is
 * truncated. Assuming the remaining sequence is not too short, the record is written to file by the
 * {@link Truncator} class that controls these objects via the {@link FastqPairParser} class.
 */
class PotentiallyTruncatedFastQRecord {
    private final String name;
    private  String sequence;
    private  String quality;

    private static String ligationSequence = null;
    private static String restrictionSequence = null;
    private static String danglingSequence = null;

    private boolean maybeDangling = false;


    static void setLigationSequence(String seq) {
        ligationSequence = seq;
    }

    static void setRestrictionSequence(String seq) {
        restrictionSequence = seq;
    }

    static void setDanglingSequence(String seq) { danglingSequence = seq; }

    public boolean isMaybeDangling() {
        return this.maybeDangling;
    }

    int getLen() { return sequence.length(); }



    PotentiallyTruncatedFastQRecord(FastqRecord fqr) {
        this.name = fqr.getReadName();
        this.sequence = fqr.getReadString();
        this.quality = fqr.getBaseQualityString();
    }

    /**
     * truncates sequence and quality lines if the {@link #ligationSequence} is found, and returns true, otherwise
     * returns false
     * @return true if truncation was performed.
     */
    boolean truncateIfLigationSiteFound() {
        int j = sequence.indexOf(PotentiallyTruncatedFastQRecord.danglingSequence);
        if (j == 0) {
            this.maybeDangling = true;
        }
        int i = sequence.indexOf(PotentiallyTruncatedFastQRecord.ligationSequence);
        if (i < 0) {
            return false;  // we did not find the ligation sequence
        }
        sequence = sequence.substring(0, i) + PotentiallyTruncatedFastQRecord.restrictionSequence;
        int len = sequence.length();
        quality = quality.substring(0, len);
        return true;
    }

    /** Write the current record to the file handle specified by out. */
    void writeToStream(BufferedWriter out) throws IOException {
        out.write("@" + name + "\n");
        out.write(sequence + "\n");
        out.write("+\n");
        out.write(quality + "\n");
    }


}
