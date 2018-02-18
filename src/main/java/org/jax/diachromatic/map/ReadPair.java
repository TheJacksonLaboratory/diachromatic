package org.jax.diachromatic.map;


import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashSet;
import java.util.Set;

import static org.jax.diachromatic.map.QCCode.*;

/**
 * This class represents a pair: one forward and one reverse read.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-06)
 */
public class ReadPair {
    private static final Logger logger = LogManager.getLogger();
    /** First (forward) read in a read pair.*/
    private final SAMRecord forwardRead;
    /** Second (reverse) read in a read pair.*/
    private final SAMRecord reverseRead;
    /** A set of Q/C criteria that this read pair did NOT pass. */
    private Set<QCCode> QCcodes;
    /** Largest allowable size of the insert of a read pair.*/
    static int UPPER_SIZE_THRESHOLD = 1500;
    /** Smallest allowable size of the insert of a read pair.*/
    static int LOWER_SIZE_THRESHOLD = 100;
    /**Length threshold in nucleotides for the end of a read being near to a restriction fragment/ligation sequence*/
    static  int DANGLING_THRESHOLD = 7;
    /** Tag to use to mark invalid reads to output to BAM file. */
    private final static String BADREAD_ATTRIBUTE="YY";
    /** Tag to mark self ligation/circularization. */
    private final static String SELF_LIGATION_TAG="SL";
    /** Tag to mark dangling end. */
    private final static String DANGLING_END_TAG="DE";
    /** Tag to same fragment internal reads. */
    private final static String SAME_INTERNAL_TAG="SI";
    /** Tag religation reads. */
    private final static String RELIGATION_TAG="RL";
    /** Tag contiguous reads. */
    private final static String CONTIGUOUS_TAG="CT";
    /** Tag for reads with too high size. */
    private final static String INSERT_TOO_BIG_TAG="TB";
    /** Tag for reads with too small size. */
    private final static String INSERT_TOO_SMALL_TAG="TS";


    ReadPair(SAMRecord f, SAMRecord r) {
        forwardRead = f;
        reverseRead = r;
        QCcodes = new HashSet<>();
    }

    // static methods to adjust threshold
    public static void setUpperSizeThreshold(int threshold) {
        UPPER_SIZE_THRESHOLD = threshold;
    }

    public static void setLowerSizeThreshold(int threshold) {
        LOWER_SIZE_THRESHOLD = threshold;
    }


    SAMRecord forward() {
        return forwardRead;
    }

    SAMRecord reverse() {
        return reverseRead;
    }

    Set<QCCode> getErrorCodes() {
        return QCcodes;
    }

    boolean isUnmapped() {
        return QCcodes.contains(READPAIR_UNMAPPED);
    }

    boolean isMultimapped() {
        return QCcodes.contains(READPAIR_MULTIMAPPED);
    }


    /**
     * This function checks if the two reads are on the sam chromosome; if not, it
     * sets the {@code CT} user-defined attribute of the reads to {@code TRANS}.
     * @param digestPair The pair of digests corresponding to the read pair. If the
     * reads are on the same chromosome, it decides whether they are {@code CLOSE}
     * or {@code FAR} and sets the {@code CT} tag accordingly.
     */
    public void characterizeReadSeparation(DigestPair digestPair) {
        if (!forwardRead.getReferenceName().equals(reverseRead.getReferenceName())) {
            // identify ditags on different chromosomes
            forwardRead.setAttribute("CT", "TRANS");
            reverseRead.setAttribute("CT", "TRANS");
            return;
        }
        // maximum possible insert size is used for determining distance of separation between fragments
        int max_possible_insert_size = digestPair.getMaximumPossibleInsertSize();
        // calculate the effective size of the insert depending on whether read 1 is mapped upstream of read 2 or vice versa
        int effective_size = forwardRead.getAlignmentStart() < reverseRead.getAlignmentStart() ?
                digestPair.reverse().getEndpos() - digestPair.forward().getStartpos() - max_possible_insert_size :
                digestPair.forward().getEndpos() - digestPair.reverse().getStartpos() - max_possible_insert_size;
        // decide whether the reads are close or far.
        if (effective_size > 10_000) {
            forwardRead.setAttribute("CT", "FAR");
            reverseRead.setAttribute("CT", "FAR");
        } else {
            forwardRead.setAttribute("CT", "CLOSE");
            reverseRead.setAttribute("CT", "CLOSE");
        }
    }

    /**
     *  <From: Wingett S et al. HiCUP: pipeline for mapping and processing Hi-C data. F1000Research 2015, 4:1310>
     *      The Hi-C protocol does not prevent entirely two adjacent restriction fragments re-ligating,
     *  but HiCUP discards such di-tags since they provide no useful three-dimensional proximity information.
     *  Similarly, multiple fragments could re-ligate forming a contig, but here paired reads will not map to
     *  adjacent genomic restriction fragments
     * This function is called if the two reads are on different fragments that are not direct neighbors. If they are located
     * within one expected fragment size, then they are contiguous sequences that were not properly digested.
     * The test demands that the contig size is above the lower threshold and below the upper threshold.
     * @return
     */
    boolean contiguous() {
        if (! forwardRead.getReferenceName().equals(reverseRead.getReferenceName())) {
            return false; // reads not on same chromosome, therefore, not contiguous
        }
//        logger.trace(String.format("contiguosus check. read1 is %s:%d-%d",readF.getReferenceName(),readF.getAlignmentStart(),readF.getAlignmentEnd()));
//        logger.trace(String.format("contiguosus check. read2 is %s:%d-%d",readR.getReferenceName(),readR.getAlignmentStart(),readR.getAlignmentEnd()));
//        logger.trace("LOWER_SIZE_THRESHOLD="+LOWER_SIZE_THRESHOLD);
//        logger.trace("readR.getAlignmentEnd() - readF.getAlignmentStart()="+(readR.getAlignmentEnd() - readF.getAlignmentStart()));
//        logger.trace("readF.getAlignmentEnd() - readR.getAlignmentStart()="+(readF.getAlignmentEnd() - readR.getAlignmentStart()));
        int contigsize=Math.max(reverseRead.getAlignmentStart() - forwardRead.getAlignmentStart(),
                forwardRead.getAlignmentStart() - reverseRead.getAlignmentStart());
        if  (contigsize >  LOWER_SIZE_THRESHOLD && contigsize < UPPER_SIZE_THRESHOLD) {
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, CONTIGUOUS_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, CONTIGUOUS_TAG);
            QCcodes.add(CONTIGUOUS);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Analyzes whether both reads are on the same restriction fragment (i.e., {@link Digest} object). If so,
     * it also analyzes which type of artifact is present -- dangling, self-ligation, or same internal.
     * @param digestPair The digest pair that corresponds to this read pair
     * @return true if both reads are located on the same digest
     */
    public boolean bothReadsLocatedOnSameRestrictionFragment(DigestPair digestPair) {
        if (!digestPair.forward().equals(digestPair.reverse())) {
            return false;
        }
        // if we get here, we know both reads are on the same restriction fragment.
        // dangling can apply to both internal and same dangling ends
        boolean dangling = danglingEnd(digestPair.forward());
        boolean selfligated = selfLigation();

        if (dangling && selfligated) {
            QCcodes.add(CIRCULARIZATION_DANGLING);
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, RELIGATION_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, RELIGATION_TAG);
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, DANGLING_END_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, DANGLING_END_TAG);
        } else if (dangling) {
            QCcodes.add(SAME_DANGLING_END);
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, DANGLING_END_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, DANGLING_END_TAG);
        } else if (selfligated) {
            QCcodes.add(CIRCULARIZATION_INTERNAL);
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, RELIGATION_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, RELIGATION_TAG);
        } else {
            // if we get here, we have reads from the same digest that are not
            // circularized and are not dangling end, so
            // they must be same_internal
            forwardRead.setAttribute(BADREAD_ATTRIBUTE, SAME_INTERNAL_TAG);
            reverseRead.setAttribute(BADREAD_ATTRIBUTE, SAME_INTERNAL_TAG);
            QCcodes.add(SAME_INTERNAL);
        }
        return true;
    }



    /**
     * If ditags are on the same restriction fragment (which MUST be checked before calling this
     * function), but not circularized and if the mapped end of one of the reads is near to the
     * end of a restriction fragment, this is termed a dangling end. Note that we only need to
     * check one Digest since by definition the reads have been found to both map to the same
     * fragment.
     * @param digest The digest that corresponds to the two reads of this readpair (same for both).
     * @return
     */
    private boolean danglingEnd(Digest digest) {
        // if we get here we know the digests are the same; take the forward one arbitrarily
        int startpos=digest.getStartpos();
        int endpos=digest.getEndpos();

        return  ( Math.abs(forwardRead.getAlignmentStart() - startpos) < DANGLING_THRESHOLD ||
                Math.abs(forwardRead.getAlignmentEnd() - endpos) < DANGLING_THRESHOLD ||
                Math.abs(reverseRead.getAlignmentStart() - startpos) < DANGLING_THRESHOLD ||
                Math.abs(reverseRead.getAlignmentEnd() - endpos) < DANGLING_THRESHOLD);
    }


    /**
     *  Adjacent fragments have the same orientation and thus the reads have opposite orientation
     *  We know the fragments are adjacent because their fragment numbers differ by 1. Recall that
     *  the {@link org.jax.diachromatic.command.DigestCommand} assigns fragments a fragment number so that
     *  adjacent fragments are number i and i+1.
     * @param digestPair
     * @return
     */
    boolean religation(DigestPair digestPair) {
        return  (Math.abs(digestPair.reverse().getFragmentNumber() - digestPair.forward().getFragmentNumber()) == 1)  &&
                (forwardRead.getReadNegativeStrandFlag() != reverseRead.getReadNegativeStrandFlag());
    }

    /** This function gets called if we cannot find valid digests for this readpair. */
    void setInvalidDigest() {
        QCcodes.add(COULD_NOT_ASSIGN_TO_DIGEST);
    }



    /**
     * Check if a fragment self ligates (circularizes). The sequence insert spans the
     * ligation site. Mapping the reads "flips" them, so that read 1 is before read2 and points in
     * the opposite direction. Vice versa if read2 is before read 1.
     * @return true if this read pair shows self-ligation
     */
    private boolean selfLigation() {
        // The following if statement is true if the reads are circularized
        return ( (forward().getAlignmentStart() < reverse().getAlignmentStart() &&
                forward().getReadNegativeStrandFlag() &&
                (!reverse().getReadNegativeStrandFlag()) )
                ||
                (reverse().getAlignmentStart() < forward().getAlignmentStart() &&
                        !forward().getReadNegativeStrandFlag() &&
                        reverse().getReadNegativeStrandFlag()) );
    }


    public boolean hasValidInsertSize(DigestPair digestPair) {
        int insertSize = getCalculatedInsertSize(digestPair);
        if (insertSize > UPPER_SIZE_THRESHOLD) {
            QCcodes.add(INSERT_TOO_LONG);
            logger.trace(String.format("Read has insert size of %d which is above the upper trheshold of %d nt", insertSize, UPPER_SIZE_THRESHOLD));
            forward().setAttribute(BADREAD_ATTRIBUTE, INSERT_TOO_BIG_TAG);
            reverse().setAttribute(BADREAD_ATTRIBUTE, INSERT_TOO_BIG_TAG);
            return false;
        } else if (insertSize < LOWER_SIZE_THRESHOLD) {
            QCcodes.add(INSERT_TOO_SHORT);
            forward().setAttribute(BADREAD_ATTRIBUTE, INSERT_TOO_SMALL_TAG);
            reverse().setAttribute(BADREAD_ATTRIBUTE, INSERT_TOO_SMALL_TAG);
            return false;
        } else {
            return true;
        }
    }



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
            logger.trace("calculated insert size: " + (end-sta+1));
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
     * add the corresponding enumeration constant from {@link QCCode} and return false. There are two
     * things that can go wrong -- either one or both reads could not be mapped, or one or both reads were mapped
     * to more than one locus in the genome. The XS attribute is used in bowtie2 to indicate that a read has
     * been multimappd
     *
     * @return true if both reads could be uniquely mapped.
     */
    boolean readPairUniquelyMapped() {
        if (forward().getReadUnmappedFlag()) {
            //read 1 could not be aligned
            QCcodes.add(READ_1_UNMAPPED);
            QCcodes.add(READPAIR_UNMAPPED);
            if (reverse().getReadUnmappedFlag()) {
                QCcodes.add(READ_2_UNMAPPED);
            }
            return false;
        } else if (reverse().getReadUnmappedFlag()) {
           // note read1 must be OK if we get here...
            QCcodes.add(READ_2_UNMAPPED);
            QCcodes.add(READPAIR_UNMAPPED);
            return false;
        } else if (forward().getAttribute("XS") != null) {
            // Now look for multimapped reads.
            // If a read has an XS attribute, then bowtie2 multi-mapped it.
            QCcodes.add(READ_1_MULTIMAPPED);
            QCcodes.add(READPAIR_MULTIMAPPED);
            if (reverse().getAttribute("XS") != null) {
                QCcodes.add(READ_2_MULTIMAPPED);
            }
            QCcodes.add(READPAIR_MULTIMAPPED);
            return false;
        } else if (reverse().getAttribute("XS") != null) {
            // note if we are here, read1 was not multimapped
            QCcodes.add(READ_2_MULTIMAPPED);
            QCcodes.add(READPAIR_MULTIMAPPED);
            return false;
        }
        return true;
    }



}
