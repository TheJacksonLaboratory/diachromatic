package org.jax.diachromatic.map;


import htsjdk.samtools.SAMRecord;

import java.util.HashSet;
import java.util.Set;

import static org.jax.diachromatic.map.ErrorCode.*;

/**
 * This class represents a pair: one forward and one reverse read.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-06)
 */
public class ReadPair {
    /** First (forward) read in a read pair. */
    private final SAMRecord forwardRead;
    /** Second (reverse) read in a read pair. */
    private final SAMRecord reverseRead;

    private Set<ErrorCode> errorcodes;

    ReadPair(SAMRecord f, SAMRecord r) {
        forwardRead=f;
        reverseRead=r;
        errorcodes=new HashSet<>();
    }


    SAMRecord forward() {
        return forwardRead;
    }


    SAMRecord reverse() {
        return reverseRead;
    }

    Set<ErrorCode> getErrorCodes(){ return errorcodes; }

    boolean isUnmapped() { return errorcodes.contains(READPAIR_UNMAPPED);}
    boolean isMultimapped(){ return errorcodes.contains(READPAIR_MULTIMAPPED);}




    /**
     * Mapped reads always "point towards" the ligation sequence. We can infer that the actualy (physical) size of the
     * insert goes from the 5' end of a read to the ligation sequence (for each read of the ditag). We calculate this
     * size and will filter out reads whose size is substantially above what we expect given the reported experimental
     * size selection step
     *
     * @param digestPair  The digest pair that corresponds to this readpair.
     * @return the insert size of chimeric read.
     */
    int getCalculatedInsertSize(DigestPair digestPair) {
        SAMRecord readF = forward();
        SAMRecord readR = reverse();
        if(!digestPair.forward().equals(digestPair.reverse())) {
            int distF, distR;
            if (readF.getReadNegativeStrandFlag()) { // readF is on the negative strand
                distF = readF.getAlignmentEnd() - digestPair.forward().getStartpos() + 1;
            } else {
                distF = digestPair.forward().getEndpos() - readF.getAlignmentStart() + 1;
            }
            if (readR.getReadNegativeStrandFlag()) { // readR is on the negative strand
                distR = readR.getAlignmentEnd() - digestPair.reverse().getStartpos() + 1;
            } else {
                distR = digestPair.reverse().getEndpos() - readR.getAlignmentStart() + 1;
            }
            return distF + distR;
        } else { // if both reads map to the same restriction fragment
            int sta=Math.min(Math.min(readF.getAlignmentStart(),readF.getAlignmentEnd()),Math.min(readR.getAlignmentStart(),readR.getAlignmentEnd()));
            int end=Math.max(Math.max(readF.getAlignmentStart(),readF.getAlignmentEnd()),Math.max(readR.getAlignmentStart(),readR.getAlignmentEnd()));
            System.out.println("XXX: " + (end-sta+1));
            return end-sta+1;
        }
    }



    /**
     * This function is called if the forward and reverse reads were found to be a valid pair
     * by {@link #readPairUniquelyMapped()}. The function adjusts the SAM flags of each read to
     * indicate that they are a valid read pair. Note that client code must call this algorithm after
     * determining that the reads should be paired, it is not done automatically by
     * {@link #readPairUniquelyMapped()}.
     */
    void pairReads() {
        // This read pair is valid
        // We therefore need to add corresponding bits to the SAM flag
        forward().setFirstOfPairFlag(true);
        reverse().setSecondOfPairFlag(true);
        // Now set the flag to indicate it is paired end data
        forward().setReadPairedFlag(true);// 0x1
        forward().setProperPairFlag(true);//0x2
        reverse().setReadPairedFlag(true);
        reverse().setProperPairFlag(true);
        // Indicate if inputSAMfiles is on the reverse strand
        forward().setMateNegativeStrandFlag(reverse().getReadNegativeStrandFlag());
        reverse().setMateNegativeStrandFlag(forward().getReadNegativeStrandFlag());

        // Set which reads are which in the inputSAMfiles
        forward().setFirstOfPairFlag(true);
        reverse().setSecondOfPairFlag(true);
        // Set the RNEXT and PNEXT values
        // If the reference indices are the same, then the following should print "="
        forward().setMateReferenceIndex(reverse().getReferenceIndex());
        reverse().setMateReferenceIndex(forward().getReferenceIndex());
        forward().setMateAlignmentStart(reverse().getAlignmentStart());
        reverse().setMateAlignmentStart(forward().getAlignmentStart());
    }



    /**
     * Determine if both reads from a paired-end could be uniquely mapped. If so, return true. If not,
     * add the corresponding enumeration constant from {@link ErrorCode} and return false. There are two
     * things that can go wrong -- either one or both reads could not be mapped, or one or both reads were mapped
     * to more than one locus in the genome. The XS attribute is used in bowtie2 to indicate that a read has
     * been multimappd
     *
     * @return true if both reads could be uniquely mapped.
     */
    boolean readPairUniquelyMapped() {
        if (forward().getReadUnmappedFlag()) {
            //read 1 could not be aligned
            errorcodes.add(READ_1_UNMAPPED);
            errorcodes.add(READPAIR_UNMAPPED);
            if (reverse().getReadUnmappedFlag()) {
                errorcodes.add(READ_2_UNMAPPED);
            }
            return false;
        } else if (reverse().getReadUnmappedFlag()) {
           // note read1 must be OK if we get here...
            errorcodes.add(READ_2_UNMAPPED);
            errorcodes.add(READPAIR_UNMAPPED);
            return false;
        } else if (forward().getAttribute("XS") != null) {
            // Now look for multimapped reads.
            // If a read has an XS attribute, then bowtie2 multi-mapped it.
            errorcodes.add(READ_1_MULTIMAPPED);
            errorcodes.add(READPAIR_MULTIMAPPED);
            if (reverse().getAttribute("XS") != null) {
                errorcodes.add(READ_2_MULTIMAPPED);
            }
            errorcodes.add(READPAIR_MULTIMAPPED);
            return false;
        } else if (reverse().getAttribute("XS") != null) {
            // note if we are here, read1 was not multimapped
            errorcodes.add(READ_2_MULTIMAPPED);
            errorcodes.add(READPAIR_MULTIMAPPED);
            return false;
        }
        return true;
    }



}
