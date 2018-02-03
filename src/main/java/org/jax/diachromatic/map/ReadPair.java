package org.jax.diachromatic.map;


import htsjdk.samtools.SAMRecord;

import java.util.HashSet;
import java.util.Set;

import static org.jax.diachromatic.map.ErrorCode.*;

/**
 * This class represents a pair: one forward and one reverse read.
 */
public class ReadPair {
    /** First (forward) read in a read pair. */
    private final SAMRecord forwardRead;
    /** Second (reverse) read in a read pair. */
    private final SAMRecord reverseRead;

    private Set<ErrorCode> errorcodes;

    public ReadPair(SAMRecord f, SAMRecord r) {
        forwardRead=f;
        reverseRead=r;
        errorcodes=new HashSet<>();
    }


    public SAMRecord forward() {
        return forwardRead;
    }


    public SAMRecord reverse() {
        return reverseRead;
    }

    public Set<ErrorCode> getErrorCodes(){ return errorcodes; }




    /**
     * This function is called if the forward and reverse reads were found to be a valid pair
     * by {@link #readPairUniquelyMapped()}. The function adjusts the SAM flags of each read to
     * indicate that they are a valid read pair. Note that client code must call this algorithm after
     * determining that the reads should be paired, it is not done automatically by
     * {@link #readPairUniquelyMapped()}.
     */
    public void pairReads() {
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
    public boolean readPairUniquelyMapped() {
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
