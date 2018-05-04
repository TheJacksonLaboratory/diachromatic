package org.jax.diachromatic.align;


import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.jax.diachromatic.align.ErrorCode.*;

/**
 * This class represents a pair of reads R1 and R2 of an illumina paired-end run.
 * <p>
 * For Hi-C,  a read pair is classified as valid (i.e. it has emerged from a DNA-DNA interaction),
 * if it fulfills the following criteria:
 * <ul>
 * <li> Both reads were mapped uniquely to the genome. </li>
 * <li> The pair does not show characteristics of known Hi-C artifacts including: </li>
 * <ul>
 * <li> <b>Same internal</b> - Both reads were mapped to the same restriction fragment. The 3' ends of the reads face one another and no read overlaps a cutting site. </li>
 * <li> <b>Same dangling end</b> - Both reads were mapped to the same restriction fragment. The 3' ends of the reads face one another and at least one read overlaps a cutting site. </li>
 * <li> <b>Same circularized alias self-ligation</b>  - Both reads were mapped to the same restriction fragment and the 5' ends of the reads face one another. This category can be further subdivided into </li>
 * <ul>
 * <li> <b>Same circularized internal</b> - No read overlaps a cutting site.
 * <li> <b>Same circularized dangling</b> - At least read overlaps a cutting site.
 * </ul>
 * <li> <b>Re-ligation</b> - Reads align to adjacent restriction fragments. </li>
 * <li> <b>Contiguous</b> - Reads align to different fragments that are not direct neighbors. But the distance between the reads is within one expected fragment size. </li>
 * <li> <b>Wrong size</b> - The calculated length of the fragment (di-tag length) is not within the limits set by the size-selection step in the experimental protocol. </li>
 * <ul>
 * <li> <b>Too small</b> - Di-tag length is smaller than a given lower threshold. </li>
 * <li> <b>Too big</b> - Di-tag length is bigger than a given lower threshold. </li>
 * </ul>
 * </ul>
 * </ul>
 * <p>
 * The insert size or di-tag length must be calculated differently for artifacts and valid pairs.
 * For artifacts the insert size is the distance between the outermost ends of mapped reads,
 * whereas for valid pairs it is the sum of the two distances from the 5' ends of mapped reads to the next occurrence of a cutting site.
 * <p>
 * This class provides all required fields and functions to assign each pair of reads (R1,R2) unambiguously to one of the categories,
 * given a align of restriction fragments ({@link Digest}) for each and the two SAMRecords of the pair.
 * <p>
 * The categorization is done upon initialization.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-06)
 */
public class ReadPair {
    private static final Logger logger = LogManager.getLogger();
    /**
     * A set of Q/C criteria that this read pair did NOT pass.
     */
    private Set<ErrorCode> errorcodes;
    /**
     * First (forward) read in a read pair.
     */
    private final SAMRecord R1;
    /**
     * Second (reverse) read in a read pair.
     */
    private final SAMRecord R2;
    /**
     * Smallest allowable size of the insert of a read pair.
     */
    private static int LOWER_SIZE_THRESHOLD = 150;
    /**
     * Largest allowable size of the insert of a read pair.
     */
    private static int UPPER_SIZE_THRESHOLD = 800;
    /**
     * Length threshold in nucleotides for the end of a read being near to a restriction fragment/ligation sequence
     */
    private static int DANGLING_THRESHOLD = 7;
    /**
     * Tag to use to mark invalid reads to output to BAM file.
     */
    private final static String BADREAD_ATTRIBUTE = "YY";
    /**
     * Set false, if read pair is an artifact.
     */
    private boolean isValid = false;
    /**
     * A read pair belongs to one of the following categories: SL, DE, CD, CI, SI, RL, TS, TL, VP.
     */
    private String categoryTag = "NA";
    /**
     * True if read 1 is unmapped.
     */
    private boolean unmapped_read1;
    /**
     * True if read 2 is unmapped.
     */
    private boolean unmapped_read2;

    /**
     * @return True if read 1 is unmapped.
     */
    boolean isUnMappedR1() {
        return unmapped_read1;
    }

    /**
     * @return True if read 2 is unmapped.
     */
    boolean isUnMappedR2() {
        return unmapped_read2;
    }

    /**
     * True if read 1 is multimapped.
     */
    private boolean multimapped_read1;
    /**
     * True if read 2 is multimapped.
     */
    private boolean multimapped_read2;

    boolean isMultiMappedR1() {
        return multimapped_read1;
    }

    boolean isMultiMappedR2() {
        return multimapped_read2;
    }

    private boolean isPaired = true;

    boolean isPaired() {
        return this.isPaired;
    }

    /**
     * Each read pair maps either to one or two fragments, here referred to as {@link Digest}, because we are using
     * an in silico digest. A {@link DigestPair} consists of two {@link Digest} objects, one for each read. Note that
     * if both reads were mapped to the same fragment, the {@link Digest} objects are identical to each other.
     */
    private final DigestPair digestPair;

    /**
     * If both reads of the pair were mapped uniquely,
     * the pair must belong to one of the categories.
     * <p>
     * Use as follows: ReadPairCategory.SELF_LIGATION.getTag() returns "SL"
     */
    private enum ReadPairCategory {
        SAME_INTERNAL("SI"),
        DANGLING_END("DE"),
        CIRULARIZED_DANGLING("CD"),
        CIRULARIZED_INTERNAL("CI"),
        RE_LIGATION("RL"),
        CONTIGUOUS("CT"),
        INSERT_TOO_SMALL("TS"),
        INSERT_TOO_LONG("TL"),
        VALID_PAIR("VP"),
        NOT_AVAILABLE("NA");

        private String tag;

        ReadPairCategory(String readPairTypeTag) {
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

    /**
     * Constructor
     *
     * @param f         forward read
     * @param r         reverse read
     * @param digestmap a align of all digests
     * @throws DiachromaticException
     */
    ReadPair(SAMRecord f, SAMRecord r, Map<String, List<Digest>> digestmap) throws DiachromaticException {

        R1 = f;
        R2 = r;

        errorcodes = new HashSet<>();
        this.digestmap = digestmap;

        // check if both reads could be mapped
        unmapped_read1 = false;
        unmapped_read2 = false;
        multimapped_read1 = false;
        multimapped_read2 = false;
        if (R1.getReadUnmappedFlag()) {
            //logger.trace("Unmapped!");
            unmapped_read1 = true;
            this.isPaired = false;
        }
        if (R2.getReadUnmappedFlag()) {
            unmapped_read2 = true;
            this.isPaired = false;
        }
        // check if both reads could be uniquely mapped
        if (R1.getAttribute("XS") != null) {
            multimapped_read1 = true;
            this.isPaired = false;
        }
        if (R2.getAttribute("XS") != null) {
            multimapped_read2 = true;
            this.isPaired = false;
        }
        // check if both reads are not on random chromosomes
        if (R1.getReferenceName().contains("_") || R2.getReferenceName().contains("_")) {
            this.isPaired = false;
        }


        if (this.isPaired) {

            // pair reads, if both reads could be mapped uniquely
            this.pairReads();

            // try to find restriction digests that match the read pair
            this.digestPair = getDigestPair(this);
            if (this.digestPair == null) {
                logger.trace("invalidDigest");
                this.setInvalidDigest();
            }

            // categorize ReadPair
            this.categorizeReadPair();
            if (!this.isValid) {
                this.R1.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
                this.R2.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
            }

            // calculate insert sizes (di-tag length)
//            if (this.isValid) {
//                logger.trace("Insert size of valid pair:" + this.getCalculatedInsertSize(digestPair) + "\t" + R1.getReadNegativeStrandFlag() + "\t" + R2.getReadNegativeStrandFlag());
//            }
        } else {
            this.digestPair=null;
        }
    }

    /**
     * Mark this read pair as valid with:
     */
    private void setValid() {
        this.isValid = true;
    }

    /**
     * Check if read pair as valid with:
     */
    public boolean isValid() {
        return this.isValid;
    }

    private void setCategoryTag(String categoryTag) {
        this.categoryTag = categoryTag;
    }

    String getCategoryTag() {
        return this.categoryTag;
    }

    // static methods to adjust threshold
    public static void setUpperSizeThreshold(int threshold) {
        UPPER_SIZE_THRESHOLD = threshold;
    }

    public static void setLowerSizeThreshold(int threshold) {
        LOWER_SIZE_THRESHOLD = threshold;
    }


    SAMRecord forward() {
        return R1;
    }

    SAMRecord reverse() {
        return R2;
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
     *
     * @param digestPair The pair of digests corresponding to the read pair. If the
     *                   reads are on the same chromosome, it decides whether they are {@code CLOSE}
     *                   or {@code FAR} and sets the {@code CT} tag accordingly.
     */
    public void characterizeReadSeparation(DigestPair digestPair) {
        if (!R1.getReferenceName().equals(R2.getReferenceName())) {
            // identify ditags on different chromosomes
            R1.setAttribute("CT", "TRANS");
            R2.setAttribute("CT", "TRANS");
            return;
        }
        // maximum possible insert size is used for determining distance of separation between fragments
        int max_possible_insert_size = digestPair.getMaximumPossibleInsertSize();
        // calculate the effective size of the insert depending on whether read 1 is mapped upstream of read 2 or vice versa
        int effective_size = R1.getAlignmentStart() < R2.getAlignmentStart() ?
                digestPair.reverse().getEndpos() - digestPair.forward().getStartpos() - max_possible_insert_size :
                digestPair.forward().getEndpos() - digestPair.reverse().getStartpos() - max_possible_insert_size;
        // decide whether the reads are close or far.
        if (effective_size > 10_000) {
            R1.setAttribute("CT", "FAR");
            R2.setAttribute("CT", "FAR");
        } else {
            R1.setAttribute("CT", "CLOSE");
            R2.setAttribute("CT", "CLOSE");
        }
    }

    /**
     * <From: Wingett S et al. HiCUP: pipeline for mapping and processing Hi-C data. F1000Research 2015, 4:1310>
     * The Hi-C protocol does not prevent entirely two adjacent restriction fragments re-ligating,
     * but HiCUP discards such di-tags since they provide no useful three-dimensional proximity information.
     * Similarly, multiple fragments could re-ligate forming a contig, but here paired reads will not align to
     * adjacent genomic restriction fragments
     * This function is called if the two reads are on different fragments that are not direct neighbors. If they are located
     * within one expected fragment size, then they are contiguous sequences that were not properly digested.
     * The test demands that the contig size is above the lower threshold and below the upper threshold.
     *
     * @return true of the two reads are contiguous (artefact!)
     */
    boolean isContiguous() {
        if (!R1.getReferenceName().equals(R2.getReferenceName())) {
            return false; // reads not on same chromosome, therefore, not contiguous
        }
        int contigsize = Math.max(R2.getAlignmentStart() - R1.getAlignmentStart(),
                R1.getAlignmentStart() - R2.getAlignmentStart());
        if (contigsize > LOWER_SIZE_THRESHOLD && contigsize < UPPER_SIZE_THRESHOLD) {  // TODO: Does this make sense?
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
     * check one Digest since by definition the reads have been found to both align to the same
     * fragment.
     *
     * @return true if the two reads of on the same restriction frag with a dangling end (artefact!)
     */
    boolean danglingEnd() {
        if (Math.abs(R1.getAlignmentStart() - digestPair.forward().getStartpos()) < DANGLING_THRESHOLD ||
                Math.abs(R1.getAlignmentStart() - digestPair.forward().getEndpos()) < DANGLING_THRESHOLD ||
                Math.abs(R2.getAlignmentStart() - digestPair.forward().getStartpos()) < DANGLING_THRESHOLD ||
                Math.abs(R2.getAlignmentStart() - digestPair.forward().getEndpos()) < DANGLING_THRESHOLD) {
            errorcodes.add(SAME_DANGLING_END);
            return true;
        } else {
            return false;
        }
    }


    /**
     * Adjacent fragments have the same orientation and thus the reads have opposite orientation
     * We know the fragments are adjacent because their fragment numbers differ by 1. Recall that
     * the {@link org.jax.diachromatic.command.DigestCommand} assigns fragments a fragment number so that
     * adjacent fragments are number i and i+1.
     *
     * @return true if the two read fragments are religated
     */
    boolean religation() {
        if ((Math.abs(digestPair.reverse().getFragmentNumber() - digestPair.forward().getFragmentNumber()) == 1) &&
                (R1.getReadNegativeStrandFlag() != R2.getReadNegativeStrandFlag())) {
            errorcodes.add(RELIGATION);
            return true;
        } else {
            return false;
        }
    }

    /**
     * This function gets called if we cannot find valid digests for this readpair.
     */
    private void setInvalidDigest() {
        errorcodes.add(COULD_NOT_ASSIGN_TO_DIGEST);
    }

    /**
     * Check if a fragment self ligates (circularizes). The sequence insert spans the
     * ligation site. Mapping the reads "flips" them, so that read 1 is before read2 and points in
     * the opposite direction. Vice versa if read2 is before read 1.
     *
     * @return true if this read pair shows self-ligation
     */
    boolean selfLigation() {
        if (!R1.getReferenceName().equals(R2.getReferenceName())) {
            return false; // reads not on same chromosome, therefore, no self-ligation
        }
        if ((forward().getAlignmentStart() < reverse().getAlignmentStart() &&
                forward().getReadNegativeStrandFlag() &&
                (!reverse().getReadNegativeStrandFlag()))
                ||
                (reverse().getAlignmentStart() < forward().getAlignmentStart() &&
                        !forward().getReadNegativeStrandFlag() &&
                        reverse().getReadNegativeStrandFlag())) {
            errorcodes.add(CIRCULARIZED_READ);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Check if insert size is too small.
     */
    private boolean hasTooSmallInsertSize() throws DiachromaticException{
        int insertSize = getCalculatedInsertSize();
        return insertSize < LOWER_SIZE_THRESHOLD;
    }

    /**
     * Check if insert size is too big.
     */
    private boolean hasTooBigInsertSize() throws DiachromaticException {
        int insertSize = getCalculatedInsertSize();
        return UPPER_SIZE_THRESHOLD < insertSize;
    }


    /**
     * Mapped reads always "point towards" the ligation sequence. We can infer that the actually (physical) size of the
     * insert goes from the 5' end of a read to the ligation sequence (for each read of the ditag). We calculate this
     * size and will filter out reads whose size is substantially above what we expect given the reported experimental
     * size selection step
     *
     * @return the insert size of chimeric read.
     */
    int getCalculatedInsertSize() throws DiachromaticException {
        SAMRecord readF = forward();
        SAMRecord readR = reverse();
        int distF, distR;
        if (digestPair==null) throw new DiachromaticException("Cannot calculate insert size with null digest align (e.g., un/multimapped read");
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
        } else { // if both reads align to the same restriction fragment
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
     * The insert size has to be calculated differently for valid pairs and other pairs.
     * This function will produce reasonable results for valid pairs only.
     *
     * @return Insert size for valid pairs.
     */
    int getInsertSizeOfValidPairs() {
        int insert_size = 0;
        int distR1, distR2;
        if (this.isValid) {
            if (R1.getReadNegativeStrandFlag()) {
                // R1 is on the negative and R2 on the positive strand
                distR1 = R1.getAlignmentEnd() - digestPair.forward().getStartpos() + 1;
                distR2 = digestPair.reverse().getEndpos() - R2.getAlignmentStart() + 1;
            } else {
                // R2 is on the negative and R1 on the positive strand
                distR1 = digestPair.forward().getEndpos() - R1.getAlignmentStart() + 1;
                distR2 = R2.getAlignmentEnd() - digestPair.reverse().getStartpos() + 1;
            }
            insert_size = distR1 + distR2;
        }
        return insert_size;
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
     * been multimapped.
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


    /* Helper functions for relative orientation of pairs */


    /**
     * Check the relative orientation of the pair: -> <-.
     *
     * @return true if the reads point to one another.
     */
    private boolean isInwardFacing() {
        return ((!this.R1.getReadNegativeStrandFlag() &&
                this.R2.getReadNegativeStrandFlag() &&
                this.R1.getAlignmentStart() < this.R2.getAlignmentEnd()) // F1R2
                ||
                (!this.R2.getReadNegativeStrandFlag() &&
                        this.R1.getReadNegativeStrandFlag() &&
                        this.R2.getAlignmentStart() < this.R1.getAlignmentEnd())); // F2R1
    }

    /**
     * @return F1F2, F2F1, R1R2, R2R1, F1R2, F2R1, R2F1, R1F2
     */
    String getRelativeOrientationTag() {

        String tag = "NA";

        if (R1.getReadNegativeStrandFlag() == R2.getReadNegativeStrandFlag()) {
            // both reads align to the same strand
            if (!R1.getReadNegativeStrandFlag()) {
                // both reads align to the forward strand
                if (R1.getAlignmentStart() <= R2.getAlignmentStart()) {
                    // R1 proceeds R2
                    tag = "F1F2";
                } else {
                    // R2 proceeds R1
                    tag = "F2F1";
                }
            } else {
                // both reads align to the reverse strand
                if (R1.getAlignmentEnd() <= R2.getAlignmentEnd()) {
                    // R1 proceeds R2
                    tag = "R1R2";
                } else {
                    // R2 proceeds R1
                    tag = "R2R1";
                }
            }
        } else {
            // reads align to different strands
            if (!R1.getReadNegativeStrandFlag()) {
                // R1 is mapped to the forward and R2 to the reverse strand
                if (R1.getAlignmentStart() <= R2.getAlignmentEnd()) {
                    // R1 proceeds R2
                    tag = "F1R2"; // innie
                } else {
                    // R2 proceeds R1
                    tag = "R2F1"; //outie
                }
            } else {
                // R1 is mapped to the reverse and R2 to the forward strand
                if (R1.getAlignmentEnd() <= R2.getAlignmentStart()) {
                    // R1 proceeds R2
                    tag = "R1F2"; // outie
                } else {
                    tag = "F2R1"; // innie
                }
            }
        }

        return tag;
    }


    /**
     * Check if at least one of the two reads overlaps the cutting site,
     * i.e. has a distance of at most DANGLING_THRESHOLD=7.
     *
     * @return true if this is the case.
     */
    private boolean readOverlapsCutSite() {
        int fragSta = this.digestPair.forward().getStartpos();
        int fragEnd = this.digestPair.forward().getEndpos();
        int fwdReadSta = this.R1.getAlignmentStart();
        int fwdReadEnd = this.R1.getAlignmentEnd();
        int revReadSta = this.R2.getAlignmentStart();
        int revReadEnd = this.R2.getAlignmentEnd();
        return (
                Math.abs(fragSta - fwdReadSta) < DANGLING_THRESHOLD ||
                        Math.abs(fragSta - fwdReadEnd) < DANGLING_THRESHOLD ||
                        Math.abs(fragSta - revReadSta) < DANGLING_THRESHOLD ||
                        Math.abs(fragSta - revReadEnd) < DANGLING_THRESHOLD ||
                        Math.abs(fragEnd - fwdReadSta) < DANGLING_THRESHOLD ||
                        Math.abs(fragEnd - fwdReadEnd) < DANGLING_THRESHOLD ||
                        Math.abs(fragEnd - revReadSta) < DANGLING_THRESHOLD ||
                        Math.abs(fragEnd - revReadEnd) < DANGLING_THRESHOLD);
    }


    private void categorizeReadPair() throws DiachromaticException {

        // 1: Determine category for read pair
        // -----------------------------------

        if (this.digestPair.forward().equals(this.digestPair.reverse())) {
            // both reads are mapped to the same fragment
            if (this.isInwardFacing()) {
                // reads point inwards
                if (this.readOverlapsCutSite()) {
                    // at least one read overlaps cutting site
                    setCategoryTag(ReadPairCategory.DANGLING_END.getTag());
                } else {
                    // no read overlaps cutting site
                    setCategoryTag(ReadPairCategory.SAME_INTERNAL.getTag());
                }
            } else {
                // reads point outwards
                if (this.readOverlapsCutSite()) {
                    // at least one read overlaps cutting site
                    setCategoryTag(ReadPairCategory.CIRULARIZED_DANGLING.getTag());
                } else {
                    // no read overlaps cutting site
                    setCategoryTag(ReadPairCategory.CIRULARIZED_INTERNAL.getTag());
                }
            }
        } else {
            // reads are mapped to different fragments
            if (this.religation()) {
                // reads are mapped to adjacent fragments
                setCategoryTag(ReadPairCategory.RE_LIGATION.getTag());
            } else {
                // reads are mapped to non adjacent fragments
                if (this.isContiguous()) {
                    // reads are located within a distance of an expected fragment size (between upper and lower threshold)
                    setCategoryTag(ReadPairCategory.CONTIGUOUS.getTag());
                } else {
                    // reads are located in a proper distance
                    if (this.hasTooSmallInsertSize()) {
                        setCategoryTag(ReadPairCategory.INSERT_TOO_SMALL.getTag());
                    } else if (this.hasTooBigInsertSize()) {
                        setCategoryTag(ReadPairCategory.INSERT_TOO_LONG.getTag());
                    } else {
                        // read pair has correct insert size
                        setCategoryTag(ReadPairCategory.VALID_PAIR.getTag());
                        this.setValid();
                    }
                }
            }
        }
    }


    /**
     * THIS FUNCTION WAS COPIED FROM THE CLASS Aligner IN ORDER TO GET THE FUNCTION isValid RUNNING.
     * IN CLASS Aligner THIS FUNCTION IS CALLED BY MULTIPLE TEST FUNCTIONS.
     * IN THIS CLASS THIS FUNCTION IS NOT YET TESTED.
     * ONE COULD MOVE THE TESTS TO THE TEST CLASS OF THIS CLASS,
     * BUT IN THE MID TERM THIS FUNCTION SHOULD BE MOVED TO THE CLASS DigestPair AND BE TESTED THERE.
     * <p>
     * Get the restriction fragments ({@link Digest} objects) to which the reads align. TODO do we need a different algorithm
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
    private DigestPair getDigestPair(ReadPair readpair) throws DiachromaticException {
        String chrom1 = readpair.forward().getReferenceName();
        int start1 = readpair.forward().getAlignmentStart();
        int end1 = readpair.forward().getAlignmentEnd();
        String chrom2 = readpair.reverse().getReferenceName();
        int start2 = readpair.reverse().getAlignmentStart();
        int end2 = readpair.reverse().getAlignmentEnd();
        final int OFFSET = 10;
        List<Digest> list = digestmap.get(chrom1);
        if (list == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not retrieve digest list for chromosome %s", chrom1));
        }
        int pos1 = start1 + OFFSET; // strand does not matter here.
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
        int pos2 = start2 + OFFSET;
        Digest d2 = list.stream().filter(digest -> (digest.getStartpos() <= pos2 && pos2 <= digest.getEndpos())).findFirst().orElse(null);
        if (d2 == null) {
            n_could_not_assign_to_digest++; // should never happen
            throw new DigestNotFoundException(String.format("Could not identify digest for read 2 at %s:%d-%d", chrom2, start2, end2));
        }

        return new DigestPair(d1, d2);

    }

    public Integer getForwardDigestStart() {return this.digestPair.forward().getStartpos();}
    public Integer getForwardDigestEnd() {return this.digestPair.forward().getEndpos();}
    public boolean forwardDigestIsActive() {return this.digestPair.forward().isActive();}

    public Integer getReverseDigestStart() {return this.digestPair.reverse().getStartpos();}
    public Integer getReverseDigestEnd() {return this.digestPair.reverse().getEndpos();}
    public boolean reverseDigestIsActive() {return this.digestPair.reverse().isActive();}



}