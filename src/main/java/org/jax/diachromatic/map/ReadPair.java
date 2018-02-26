package org.jax.diachromatic.map;


import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;

import java.awt.*;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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
    private static final Logger logger = LogManager.getLogger();
    /** First (forward) read in a read pair.*/
    private final SAMRecord forwardRead;
    /** Second (reverse) read in a read pair.*/
    private final SAMRecord reverseRead;
    /** A set of Q/C criteria that this read pair did NOT pass. */
    private Set<ErrorCode> errorcodes;
    /** Largest allowable size of the insert of a read pair.*/
    static int UPPER_SIZE_THRESHOLD = 1500;
    /** Smallest allowable size of the insert of a read pair.*/
    static int LOWER_SIZE_THRESHOLD = 100;
    /**Length threshold in nucleotides for the end of a read being near to a restriction fragment/ligation sequence*/
    static  int DANGLING_THRESHOLD = 7;
    /** Tag to use to mark invalid reads to output to BAM file. */
    private final static String BADREAD_ATTRIBUTE="YY";
    /** Set false, if read pair is an artifact. */
    private boolean isValid=false;
    /** A read pair belongs to one of the following categories: SL, DE, CD, CI, SI, RL, TS, TB, VP. */
    private String categoryTag="NA";
    /** Insert size */
    private int insertSize;

    /**
     * Each read pair maps either to one or two fragments.
     *
     * What is digestPair, if both reads were mapped to the same fragment?
     */
    private DigestPair digestPair;

    /**
     *
     * Use as follows: readPairCategory.SELF_LIGATION.getTag() returns "SL"
     *
     */
    private enum readPairCategory
    {
        DANGLING_END("DE"),
        CIRULARIZED_DANGLING("CD"),
        CIRULARIZED_INTERNAL("CI"),
        SAME_INTERNAL("SI"),
        RE_LIGATION("RL"),
        CONTIGUOUS("CT"),
        INSERT_TOO_SMALL("TS"),
        INSERT_TOO_BIG("TB"),
        VALID_PAIR("VP");

        private String tag;

        readPairCategory(String readPairTypeTag) {
            this.tag = readPairTypeTag;
        }

        public String getTag() {
            return tag;
        }
    }

    /**
     * Key: chromosome; value: a list of {@link Digest} objects on the chromosome.
     */
    private Map<String, List<Digest>> digestmap = null;

    private int n_could_not_assign_to_digest = 0;


    ReadPair(SAMRecord f, SAMRecord r, Map<String, List<Digest>> digestmap) throws DiachromaticException {
        forwardRead = f;
        reverseRead = r;
        errorcodes = new HashSet<>();
        this.digestmap=digestmap;

        this.digestPair = getDigestPair(this);
        // check whether we can find restriction digests that match the read pair
        if (this.digestPair == null) {
            this.setInvalidDigest();
        }

        // categorize ReadPair
        this.categorizeReadPair();
        if(!this.isValid) {
            this.forwardRead.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
            this.reverseRead.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
        }


        // set insert size (calculation of the insert size depends on the category)


    }

    /**
     * Mark this read pair as valid.
     */
    private void setValid() {
        this.isValid=true;
    }

    /**
     * Check if read pair as valid.
     */
    public boolean  isValid() {
        return this.isValid;
    }

    private void setCategoryTag(String categoryTag) {
        this.categoryTag=categoryTag;
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

    Set<ErrorCode> getErrorCodes() {
        return errorcodes;
    }

    boolean isUnmapped() {
        return errorcodes.contains(READPAIR_UNMAPPED);
    }

    boolean isMultimapped() {
        return errorcodes.contains(READPAIR_MULTIMAPPED);
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
    boolean isContiguous() {
        if (! forwardRead.getReferenceName().equals(reverseRead.getReferenceName())) {
            return false; // reads not on same chromosome, therefore, not contiguous
        }
        int contigsize=Math.max(reverseRead.getAlignmentStart() - forwardRead.getAlignmentStart(),
                forwardRead.getAlignmentStart() - reverseRead.getAlignmentStart());
        if  (contigsize >  LOWER_SIZE_THRESHOLD && contigsize < UPPER_SIZE_THRESHOLD) {
            errorcodes.add(CONTIGUOUS);
            return true;
        } else {
            return false;
        }
    }

    /**
     * If ditags are on the same restriction fragment (which MUST be checked before calling this
     * function), but not circularized and if the mapped end of one of the reads is near to the
     * end of a restriction fragment, this is termed a dangling end. Note that we only need to
     * check one Digest since by definition the reads have been found to both map to the same
     * fragment.
     * @param digestPair The digest pair that corresponds to the two reads of this readpair.
     * @return
     */
    boolean danglingEnd(DigestPair digestPair) {
        if  ( Math.abs(forwardRead.getAlignmentStart() - digestPair.forward().getStartpos()) < DANGLING_THRESHOLD ||
                Math.abs(forwardRead.getAlignmentStart() - digestPair.forward().getEndpos()) < DANGLING_THRESHOLD ||
                Math.abs(reverseRead.getAlignmentStart() - digestPair.forward().getStartpos()) < DANGLING_THRESHOLD ||
                Math.abs(reverseRead.getAlignmentStart() - digestPair.forward().getEndpos()) < DANGLING_THRESHOLD) {
            errorcodes.add(SAME_DANGLING_END);
            return true;
        } else {
            return false;
        }
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
        if  ( (Math.abs(digestPair.reverse().getFragmentNumber() - digestPair.forward().getFragmentNumber()) == 1)  &&
                (forwardRead.getReadNegativeStrandFlag() != reverseRead.getReadNegativeStrandFlag()) ) {
            errorcodes.add(RELIGATION);
            return true;
        } else {
            return false;
        }
    }

    /** This function gets called if we cannot find valid digests for this readpair. */
    void setInvalidDigest() {
        errorcodes.add(COULD_NOT_ASSIGN_TO_DIGEST);
    }

    /**
     * Check if a fragment self ligates (circularizes). The sequence insert spans the
     * ligation site. Mapping the reads "flips" them, so that read 1 is before read2 and points in
     * the opposite direction. Vice versa if read2 is before read 1.
     * @return true if this read pair shows self-ligation
     */
    boolean selfLigation() {
        if (! forwardRead.getReferenceName().equals(reverseRead.getReferenceName())) {
            return false; // reads not on same chromosome, therefore, no self-ligation
        }
        if ( (forward().getAlignmentStart() < reverse().getAlignmentStart() &&
                forward().getReadNegativeStrandFlag() &&
                (!reverse().getReadNegativeStrandFlag()) )
                ||
                (reverse().getAlignmentStart() < forward().getAlignmentStart() &&
                        !forward().getReadNegativeStrandFlag() &&
                        reverse().getReadNegativeStrandFlag()) ) {
            errorcodes.add(CIRCULARIZED_READ);
            return true;
        }  else {
            return false;
        }
    }

    /**
     * Check if insert size is too small.
     */
     private boolean hasTooSmallInsertSize() {
         int insertSize = getCalculatedInsertSize(this.digestPair);
         if (insertSize < LOWER_SIZE_THRESHOLD) {
             return true;
         }
         else {
             return false;
         }
     }

    /**
     * Check if insert size is too big.
     */
    private boolean hasTooBigInsertSize() {
        int insertSize = getCalculatedInsertSize(this.digestPair);
        if (UPPER_SIZE_THRESHOLD < insertSize) {
            return true;
        }
        else {
            return false;
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
        int distF, distR;
        if (!digestPair.forward().equals(digestPair.reverse())) {
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
            // Handle size calculation for circularized fragments. Needs to be calculated as above.
            // If the read mapped to the reverse strand comes before the read mapped to the forward strand, calculate insert size as above
            int insert_size;
            logger.trace("Both reads are on the same fragment.");

            if (!readF.getMateNegativeStrandFlag() && readF.getAlignmentStart() > readR.getAlignmentStart()) { // readF is on the positive strand and comes after the other read
                distF = digestPair.forward().getEndpos() - readF.getAlignmentStart() + 1;
                distR = readR.getAlignmentEnd() - digestPair.reverse().getStartpos() + 1;
                insert_size = distF + distR;
            } else if (readF.getMateNegativeStrandFlag() && readR.getAlignmentStart() > readF.getAlignmentStart()) {
                distF = readF.getAlignmentEnd() - digestPair.forward().getStartpos() + 1;
                distR = digestPair.reverse().getEndpos() - readR.getAlignmentStart() + 1;
                insert_size = distF + distR;
            } else {
                int sta = Math.min(Math.min(readF.getAlignmentStart(), readF.getAlignmentEnd()), Math.min(readR.getAlignmentStart(), readR.getAlignmentEnd()));
                int end = Math.max(Math.max(readF.getAlignmentStart(), readF.getAlignmentEnd()), Math.max(readR.getAlignmentStart(), readR.getAlignmentEnd()));
                logger.trace("calculated insert size: " + (end - sta + 1));
                insert_size = end - sta + 1;
            }
            return insert_size;
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

    /**
     * Check the relative orientation of the pair.
     *
     * @return true if the reads point to one another.
     */
    boolean isFacingPair() {
        if((!this.forwardRead.getReadNegativeStrandFlag() && this.reverseRead.getReadNegativeStrandFlag() &&  this.forwardRead.getAlignmentStart()<this.reverseRead.getAlignmentEnd())
           ||
           (!this.reverseRead.getReadNegativeStrandFlag() && this.forwardRead.getReadNegativeStrandFlag() &&  this.reverseRead.getAlignmentStart()<this.forwardRead.getAlignmentEnd()))
        {
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Check if at least one of the two reads overlaps the cutting site,
     * i.e. has a distance of at most DANGLING_THRESHOLD=7.
     *
     * @return true if this is the case.
     */
    boolean readOverlapsCutSite() {
        int fragSta=this.digestPair.forward().getStartpos();
        int fragEnd=this.digestPair.forward().getEndpos();
        int fwdReadSta=this.forwardRead.getAlignmentStart();
        int fwdReadEnd=this.forwardRead.getAlignmentEnd();
        int revReadSta=this.reverseRead.getAlignmentStart();
        int revReadEnd=this.reverseRead.getAlignmentEnd();
        if(
            Math.abs(fragSta-fwdReadSta) < DANGLING_THRESHOLD ||
            Math.abs(fragSta-fwdReadEnd) < DANGLING_THRESHOLD ||
            Math.abs(fragSta-revReadSta) < DANGLING_THRESHOLD ||
            Math.abs(fragSta-revReadEnd) < DANGLING_THRESHOLD ||
            Math.abs(fragEnd-fwdReadSta) < DANGLING_THRESHOLD ||
            Math.abs(fragEnd-fwdReadEnd) < DANGLING_THRESHOLD ||
            Math.abs(fragEnd-revReadSta) < DANGLING_THRESHOLD ||
            Math.abs(fragEnd-revReadEnd) < DANGLING_THRESHOLD)
        {
            return true;
        } else {
            return false;
        }
    }


    private void categorizeReadPair() throws DiachromaticException {

        // 1: Determine category for read pair
        // -----------------------------------

        if (this.digestPair.forward().equals(this.digestPair.reverse())) {
            // both reads are mapped to the same fragment
            if(this.isFacingPair())
            {
                // reads point inwards
                if(this.readOverlapsCutSite()) {
                    // at least one read overlaps cutting site
                    setCategoryTag(readPairCategory.DANGLING_END.getTag());
                }
                else {
                    // no read overlaps cutting site
                    setCategoryTag(readPairCategory.SAME_INTERNAL.getTag());
                }
            } else {
                // reads point outwards
                if(this.readOverlapsCutSite()) {
                    // at least one read overlaps cutting site
                    setCategoryTag(readPairCategory.CIRULARIZED_DANGLING.getTag());
                }
                else {
                    // no read overlaps cutting site
                    setCategoryTag(readPairCategory.CIRULARIZED_INTERNAL.getTag());
                }
            }
        }
        else {
            // reads are mapped to different fragments
            if(this.religation(this.digestPair)) {
                // reads are mapped to adjacent fragments
                setCategoryTag(readPairCategory.RE_LIGATION.getTag());
            }
            else {
                // reads are mapped to non adjacent fragments
                if(this.isContiguous()) {
                    // reads are located within a distance of an expected fragment size (between upper and lower threshold)
                    setCategoryTag(readPairCategory.CONTIGUOUS.getTag());
                }
                else {
                    // reads are located in a proper distance
                    if(this.hasTooSmallInsertSize()) {
                        setCategoryTag(readPairCategory.INSERT_TOO_SMALL.getTag());
                    }
                    else if(this.hasTooBigInsertSize()) {
                        setCategoryTag(readPairCategory.INSERT_TOO_BIG.getTag());
                    } else {
                        // read pair has correct insert size
                        setCategoryTag(readPairCategory.VALID_PAIR.getTag());
                        this.setValid();
                    }
                }
            }
        }
    }


    /**
     *
     * THIS FUNCTION WAS COPIED FROM THE CLASS SAMPairer IN ORDER TO GET THE FUNCTION isValid RUNNING.
     * IN CLASS SAMPairer THIS FUNCTION IS CALLED BY MULTIPLE TEST FUNCTIONS.
     * IN THIS CLASS THIS FUNCTION IS NOT YET TESTED.
     * ONE COULD MOVE THE TESTS TO THE TEST CLASS OF THIS CLASS,
     * BUT IN THE MID TERM THIS FUNCTION SHOULD BE MOVED TO THE CLASS DigestPair AND BE TESTED THERE.
     *
     * Get the restriction fragments ({@link Digest} objects) to which the reads map. TODO do we need a different algorithm
     * Note this from hicup
     * Using the terminal ends of a di-tag to position a read on the digested genome could be problematic because
     * a restiction enzyme may not cut in the same position on both strands (i.e. creates sticky ends). Filling-in
     * and truncating reads at the Hi-C junction makes this situation even more complex.
     * To overcome this problem simply the script uses the position 10 bp upstream of the start of each read when
     * assigning reads to a fragment in the digested genome.
     *
     * @param readpair Pair of reads (forward, reverse).
     * @return the corresponding {@link DigestPair} object.
     */
    DigestPair getDigestPair(ReadPair readpair) throws DiachromaticException {
        String chrom1 = readpair.forward().getReferenceName();
        int start1 = readpair.forward().getAlignmentStart();
        int end1 = readpair.forward().getAlignmentEnd();
        String chrom2 = readpair.reverse().getReferenceName();
        int start2 = readpair.reverse().getAlignmentStart();
        int end2 = readpair.reverse().getAlignmentEnd();
        return getDigestPair(chrom1, start1, end1, chrom2, start2, end2);
    }

    /**
     *
     * THIS FUNCTION WAS COPIED FROM THE CLASS SAMPairer IN ORDER TO GET THE FUNCTION isValid RUNNING.
     * IN CLASS SAMPairer THIS FUNCTION IS CALLED BY MULTIPLE TEST FUNCTIONS.
     * IN THIS CLASS THIS FUNCTION IS NOT YET TESTED.
     * ONE COULD MOVE THE TESTS TO THE TEST CLASS OF THIS CLASS,
     * BUT IN THE MID TERM THIS FUNCTION SHOULD BE MOVED TO THE CLASS DigestPair AND BE TESTED THERE.
     *
     * Using the terminal ends of a di-tag to position a read on the digested genome could be problematic because
     * a restiction enzyme may not cut in the same position on both strands (i.e. creates sticky ends). Filling-in
     * and truncating reads at the Hi-C junction makes this situation even more complex.
     * To overcome this problem we use the position 10 bp upstream of the start of each read when
     * assigning reads to a fragment in the digested genome. This approach was adopted from the HiCup script.
     *
     * @return
     */
    DigestPair getDigestPair(String chrom1, int start1, int end1, String chrom2, int start2, int end2) throws DiachromaticException {
        final int OFFSET=10;
        List<Digest> list = digestmap.get(chrom1);
        if (list == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not retrieve digest list for chromosome %s", chrom1));
        }
        int pos1=start1+OFFSET; // strand does not matter here.
        Digest d1 = list.stream().filter(digest -> (digest.getStartpos() <= pos1 && pos1 <= digest.getEndpos())).findFirst().orElse(null);
        if (d1 == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not identify digest for read 1 at %s:%d-%d", chrom1, start1, end1));
        }
        list = digestmap.get(chrom2);
        if (list == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not retrieve digest list for chromosome %s", chrom2));
        }
        int pos2=start2+OFFSET;
        Digest d2 = list.stream().filter(digest -> (digest.getStartpos() <= pos2 && pos2 <= digest.getEndpos())).findFirst().orElse(null);
        if (d2 == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not identify digest for read 2 at %s:%d-%d", chrom2, start2, end2));
        }

        return new DigestPair(d1, d2);

    }
}
