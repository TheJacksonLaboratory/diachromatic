package org.jax.diachromatic.align;

import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

/**
 * This class represents a pair of reads R1 and R2 of an Illumina paired-end run for a Hi-C library.
 * <p>
 * The two reads represent the ends of fragments that are ideally valid Hi-C hybrid fragments arising from functional
 * interactions and consisting of two pieces of DNA originating from two interacting loci. However, there are
 * various sources of artifacts that lead to un-ligated or self-ligated fragments and corresponding read pairs.
 * Furthermore, there are hybrid fragments whose size is inconsistent with parameters chosen for sonication.
 * Based on this, Diachromatic subdivides the set of uniquely mapping pairs into five categories:
 * <ul>
 * <li> <b>Un-ligated</b> - The 5' end positions of the mapped reads have distance that is consistent with fragment sizes that result from sonication. </li>
 * <li> <b>Self-ligated</b> - The calculated size of the hypothetical self-ligated fragment is smaller than a given threshold. </li>
 * <li> <b>Too short</b>  - The read pair was not categorized as un-ligated or self-ligated, but the calculated size of the corresponding hybrid fragment is too small. </li>
 * <li> <b>Too long</b> - Same as for <i>too long</i> Reads align to adjacent restriction fragments. </li>
 * <li> <b>Valid</b> - All remaining pairs are categorized as valid pairs that can be used for downstream analysis. </li>
 * </ul>
 * <p>
 * There two further features of read pairs that do not affect categorization. If one of the 5' end positions of the
 * mapped reads coincides with a restriction enzyme cutting site the read pair is referred to as <b>dangling end read
 * pairs</b>. Large numbers of dangling end pairs indicate problems re-ligation conditions. The second feature relates
 * to cross-ligation between the digested protein-DNA complexes that result in <b>trans read pairs</b> for which the two
 * reads map to different chromosomes. of a cutting site.
 * <p>
 * See documentation on read the docs or manuscript.
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
     * First (forward) and second (reverse) read in a read pair.
     */
    private SAMRecord R1=null;
    private SAMRecord R2=null;

    /**
     * These two thresholds define the range of fragment sizes that are considered to be consistent with the chosen
     * parameters for sonication. For capture Hi-C, an average sizes from 300 to 500 bp are recommended.
     *
     * Inward-pointing read pairs whose 5' end positions have a distance d smaller than the defined upper threshold will
     * be categorized as un-ligated.
     *
     * Valid read pairs for which the calculated hybrid fragment sizes d' is outside the specified range will be
     * categorized as too small or too large.
     *
     * For outward pointing read pairs, the calculated size of the hybrid fragment d' is added to d in order to calculate
     * the size of the corresponding self-ligated fragment d''. If d'' is smaller than the upper threshold, the read pair
     * will be categorized as self-ligated.
     *
     * TODO: Discuss if there should be a separate threshold for self-ligation.
     */
    private static int LOWER_SIZE_THRESHOLD = 150;
    private static int UPPER_SIZE_THRESHOLD = 800;

    /**
     * Length threshold in nucleotides for the end of a read being near to a restriction fragment/ligation sequence.
     */
    private static int DANGLING_THRESHOLD = 7;

    /**
     * Tag to use to mark invalid reads to output to BAM file.
     */
    private final static String BADREAD_ATTRIBUTE = "YY";

    /**
     * Tag to indicate relative orientation of read pair ()
     */
    private final static String ORIENTATION_ATTRIBUTE = "RO";

    /**
     * True, if read R1 or R2 is unmapped.
     */
    private boolean unmapped_R1;
    private boolean unmapped_R2;

    /**
     * @return True, if R1 is unmapped.
     */
    boolean isUnMappedR1() {
        return unmapped_R1;
    }

    /**
     * @return True, if R2 is unmapped.
     */
    boolean isUnMappedR2() {
        return unmapped_R2;
    }

    /**
     * True, if read R1 or R2 is multi-mapped.
     */
    private boolean multimapped_R1;
    private boolean multimapped_R2;

    /**
     * @return True, if R1 is multi-mapped.
     */
    boolean isMultiMappedR1() {
        return multimapped_R1;
    }

    /**
     * @return @return True, if R1 is multi-mapped.
     */
    boolean isMultiMappedR2() {
        return multimapped_R2;
    }

    /**
     * This field is set to true, if both reads of the pair can be uniquely mapped.
     */
    private boolean isPaired = true;

    /**
     * @return True, if both reads of the pair can be uniquely mapped.
     */
    boolean isPaired() {
        return this.isPaired;
    }

    /**
     * Each read pair maps either to one or two fragments, here referred to as {@link Digest}, because we are using
     * an in silico digest. A {@link DigestPair} consists of two {@link Digest} objects, one for each read. Note that
     * if both reads were mapped to the same fragment, the {@link Digest} objects are identical to each other.
     */
    private DigestPair digestPair = null;

    /**
     * A read pair belongs to one of the following categories: UL, SL, TS, TL or VP.
     */
    private String categoryTag = "NA";

    /**
     * Paired reads belong to one of the following disjoint categories.
     */
    private enum ReadPairCategory {
        UN_LIGATED("UL"),
        SELF_LIGATED("SL"),
        VALID_PAIR("VP"),
        VALID_TOO_SHORT("TS"),
        VALID_TOO_LONG("TL");

        private String tag;

        ReadPairCategory(String readPairTypeTag) {
            this.tag = readPairTypeTag;
        }

        public String getTag() {
            return tag;
        }
    }

    private void setCategoryTag(String categoryTag) {
        this.categoryTag = categoryTag;
    }

    String getCategoryTag() {
        return this.categoryTag;
    }

    /**
     * Simple constructor that is used for counting of interactions.
     *
     * @param f SAM record for R1
     * @param r SAM record for R2
     * @param digestMap Custom class for storing all digests of the genome.
     * @throws DiachromaticException Required because getDigestPair2 is used.
     */
    public ReadPair(SAMRecord f, SAMRecord r, DigestMap digestMap) throws DiachromaticException {

        R1 = f;
        R2 = r;

        // create digest pair
        this.digestPair = digestMap.getDigestPair2(f.getReferenceName(),getFivePrimeEndPosOfRead(f),r.getReferenceName(),getFivePrimeEndPosOfRead(r));
    }


    /**
     *
     * @param f SAMRecord for R1.
     * @param r SAMRecord for R2.
     * @param digestMap structure containing all digests.
     * @param lowerFragSize lower size threshold for fragments after sonication.
     * @param upperFragSize upper size threshold for fragments after sonication.
     * @param stringentUnique If true, more  stringent definition of uniquely mapped is used.
     * @throws DiachromaticException
     */
    ReadPair(SAMRecord f, SAMRecord r, DigestMap digestMap, Integer lowerFragSize, Integer upperFragSize, boolean stringentUnique) throws DiachromaticException {

        R1 = f;
        R2 = r;

        this.LOWER_SIZE_THRESHOLD = lowerFragSize;
        this.UPPER_SIZE_THRESHOLD = upperFragSize;

        // check if both reads could be mapped
        unmapped_R1 = false;
        unmapped_R2 = false;
        if (R1.getReadUnmappedFlag()) {
            unmapped_R1 = true;
            this.isPaired = false;
        }
        if (R2.getReadUnmappedFlag()) {
            unmapped_R2 = true;
            this.isPaired = false;
        }

        // check if both reads could be uniquely mapped
        multimapped_R1 = false;
        multimapped_R2 = false;

        if (R1.getAttribute("XS") != null) {
            // there is more than one alignment
            if(stringentUnique) {
                // in the stringent mode this enough to be multi-mapped
                multimapped_R1 = true;
                this.isPaired = false;
            } else {
                if(R1.getMappingQuality()<30 || (int)R1.getAttribute("AS")-(int)R1.getAttribute("XS")<10) {
                    multimapped_R1 = true;
                    this.isPaired = false;
                }
            }
        }
        if (R2.getAttribute("XS") != null) {
            // there is more than one alignment
            if(stringentUnique) {
                // in the stringent mode this enough to be multi-mapped
                multimapped_R2 = true;
                this.isPaired = false;
            } else {
                if(R2.getMappingQuality()<30 || (int)R2.getAttribute("AS")-(int)R2.getAttribute("XS")<10) {
                    multimapped_R2 = true;
                    this.isPaired = false;
                }
            }
        }

        // check if both reads are not on random chromosomes
        if (R1.getReferenceName().contains("_") || R2.getReferenceName().contains("_")) {
            this.isPaired = false;
        }

        if (this.isPaired) {

            // pair reads, if both reads could be mapped uniquely
            this.pairReads();

            // try to find restriction digests that match the read pair
            this.digestPair = digestMap.getDigestPair2(this.R1.getReferenceName(),getFivePrimeEndPosOfRead(this.R1),this.R2.getReferenceName(),getFivePrimeEndPosOfRead(this.R2));

            // categorize ReadPair
            this.categorizeReadPair();
            if (!this.getCategoryTag().equals("VP")) {
                this.R1.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
                this.R2.setAttribute(BADREAD_ATTRIBUTE, this.categoryTag);
            }

            // add attribute for relative orientation of read pair
            this.R1.setAttribute(ORIENTATION_ATTRIBUTE, this.getRelativeOrientationTag());
            this.R2.setAttribute(ORIENTATION_ATTRIBUTE, this.getRelativeOrientationTag());


        } else {
            this.digestPair=null;
        }
    }

    /**
     * @param samRecord
     * @return The genomic position that corresponds to the 5' end position of the mapped read.
     */
    private Integer getFivePrimeEndPosOfRead(SAMRecord samRecord) {
        if(!samRecord.getReadNegativeStrandFlag()) {
            return samRecord.getAlignmentStart();
        }
        else {
            return samRecord.getAlignmentEnd();
        }
    }

    public Integer getFivePrimeEndPosOfR1() {
        if(!this.R1.getReadNegativeStrandFlag()) {
            return this.R1.getAlignmentStart();
        }
        else {
            return this.R1.getAlignmentEnd();
        }
    }

    public Integer getFivePrimeEndPosOfR2() {
        if(!this.R2.getReadNegativeStrandFlag()) {
            return this.R2.getAlignmentStart();
        }
        else {
            return this.R2.getAlignmentEnd();
        }
    }

    public SAMRecord forward() {
        return R1;
    }

    public SAMRecord reverse() {
        return R2;
    }

    /**
     * This function checks if the two reads are on the sam chromosome; if not, it
     * sets the {@code CT} user-defined attribute of the reads to {@code TRANS}.
     *
     * @param digestPair The pair of digests corresponding to the read pair. If the
     *                   reads are on the same chromosome, it decides whether they are {@code CLOSE}
     *                   or {@code FAR} and sets the {@code CT} tag accordingly.
     *
     * TODO: Discuss the usefulness of the function below and discard or use.
     *
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
     * @return True, if both reads are on the same chromosome.
     */
    public boolean isTrans() {
        if(!R1.getReferenceName().equals(R2.getReferenceName())) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Mapped reads always "point towards" the ligation sequence. We can infer that the actually (physical) size of the
     * insert goes from the 5' end of a read to the ligation sequence (for each read of the ditag). We calculate this
     * size and will filter out reads whose size is substantially above what we expect given the reported experimental
     * size selection step. TODO: Revise!
     *
     * @return The insert size for chimeric hybrid fragments and also for un-ligated read pairs.
     */
    public Integer getCalculatedInsertSize() {

         SAMRecord R1 = forward();
         SAMRecord R2 = reverse();

         /*
         For un-ligated read pairs the size corresponds to the distance between the 5' end postions of the mapped reads.
           */
         if(this.getCategoryTag().equals("UL")) {
             return this.getDistanceBetweenFivePrimeEnds();
         }

         int d1;
         if(!R1.getReadNegativeStrandFlag()) {
             d1 = digestPair.forward().getEndpos() - getFivePrimeEndPosOfRead(R1) + 1;
         }
         else {
             d1 = getFivePrimeEndPosOfRead(R1) - digestPair.forward().getStartpos() + 1;
         }
         int d2;
         if(!R2.getReadNegativeStrandFlag()) {
             d2 = digestPair.reverse().getEndpos() - getFivePrimeEndPosOfRead(R2) + 1;
         }
         else {
             d2 = getFivePrimeEndPosOfRead(R2) - digestPair.reverse().getStartpos() + 1;
         }
         return d1 + d2;
     }


    /**
     * The two reads of given pairs are independently mapped to the genome as if they were single-end reads.
     * If both reads could be mapped uniquely, this function is used in order to adjust the SAM fields of the two reads
     * so as they can be recognized as mapped paired-end read.
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


    /* Helper functions for relative orientation of pairs */

    /**
     * Check the relative orientation of the pair.
     *
     * @return True, if the reads point to one another.
     */
    private boolean isInwardFacing() {
        if(getRelativeOrientationTag().equals("F1R2") || getRelativeOrientationTag().equals("F2R1")) {
            return true;
        } else {
            return false;
        }
    }


    /**
     * Check the relative orientation of the pair.
     *
     * @return True, if the reads point to opposite directions.
     */
    private boolean isOutwardFacing() {
        if(getRelativeOrientationTag().equals("R2F1") || getRelativeOrientationTag().equals("R1F2")) {
            return true;
        } else {
            return false;
        }
    }


    /**
     * Get relative orientation of the pair.
     *
     * @return F1F2, F2F1, R1R2, R2R1, F1R2, F2R1, R2F1 or R1F2.
     */
    public String getRelativeOrientationTag() {

        String tag = "NA";

        if (R1.getReadNegativeStrandFlag() == R2.getReadNegativeStrandFlag()) {
            // both reads align to the same strand
            if (!R1.getReadNegativeStrandFlag()) {
                // both reads align to the forward strand
                if (getFivePrimeEndPosOfRead(R1) <= getFivePrimeEndPosOfRead(R2)) {
                    // R1 proceeds R2
                    tag = "F1F2";
                } else {
                    // R2 proceeds R1
                    tag = "F2F1";
                }
            } else {
                // both reads align to the reverse strand
                if (getFivePrimeEndPosOfRead(R1) <= getFivePrimeEndPosOfRead(R2)) {
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
                if (getFivePrimeEndPosOfRead(R1) <= getFivePrimeEndPosOfRead(R2)) {
                    // R1 proceeds R2
                    tag = "F1R2"; // innie
                } else {
                    // R2 proceeds R1
                    tag = "R2F1"; //outie
                }
            } else {
                // R1 is mapped to the reverse and R2 to the forward strand
                if (getFivePrimeEndPosOfRead(R1) <= getFivePrimeEndPosOfRead(R2)) {
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
     * Check if at least one of the two reads overlaps the cutting site, i.e. the 5' end position has a distance of at
     * most DANGLING_THRESHOLD = 7.
     *
     * @return True, if this pair belongs to a dangling end fragment.
     */
    public boolean isDanglingEnd() {
        int fragSta = this.digestPair.forward().getStartpos();
        int fragEnd = this.digestPair.forward().getEndpos();

        int fwdReadFpep = getFivePrimeEndPosOfRead(this.R1);
        int revReadFpep = getFivePrimeEndPosOfRead(this.R2);
        return(
                Math.abs(fragSta-fwdReadFpep) <  DANGLING_THRESHOLD || Math.abs(fragEnd-fwdReadFpep) <  DANGLING_THRESHOLD ||
                Math.abs(fragSta-revReadFpep) <  DANGLING_THRESHOLD || Math.abs(fragEnd-revReadFpep) <  DANGLING_THRESHOLD);
    }


    /**
     * Assigns this read pair to one of the following disjoint categories: Self-ligated, un-ligated, wrong size or valid.
     *
     * @throws DiachromaticException
     */
    private void categorizeReadPair() throws DiachromaticException {

        if(!this.isTrans() && (this.isInwardFacing()||this.isOutwardFacing())){
            // inward and outward pointing read pairs may arise from  self- or un-ligated fragments
            if(this.isOutwardFacing()) {
                if(this.getDistanceBetweenFivePrimeEnds() < this.UPPER_SIZE_THRESHOLD) {
                    setCategoryTag(ReadPairCategory.SELF_LIGATED.getTag());
                } else {
                    if(this.hasTooSmallInsertSize()) {
                        setCategoryTag(ReadPairCategory.VALID_TOO_SHORT.getTag());
                    } else if(this.hasTooBigInsertSize()){
                        setCategoryTag(ReadPairCategory.VALID_TOO_LONG.getTag());
                    } else {
                        setCategoryTag(ReadPairCategory.VALID_PAIR.getTag()); // only if hybrid fragment has the right size it is categorized as valid
                    }
                }
            } else {
                // pair is inward facing
                if(this.getDistanceBetweenFivePrimeEnds() < this.UPPER_SIZE_THRESHOLD) {
                    setCategoryTag(ReadPairCategory.UN_LIGATED.getTag());
                } else {
                    if(this.hasTooSmallInsertSize()) {
                        setCategoryTag(ReadPairCategory.VALID_TOO_SHORT.getTag());
                    } else if(this.hasTooBigInsertSize()){
                        setCategoryTag(ReadPairCategory.VALID_TOO_LONG.getTag());
                    } else {
                        setCategoryTag(ReadPairCategory.VALID_PAIR.getTag()); // only if hybrid fragment has the right size it is categorized as valid
                    }
                }
            }
        } else {
            // trans pairs and read pairs that are pointing in the same direction cannot arise from self- or un-ligated fragments
            if(this.hasTooSmallInsertSize()) {
                setCategoryTag(ReadPairCategory.VALID_TOO_SHORT.getTag());
            } else if(this.hasTooBigInsertSize()){
                setCategoryTag(ReadPairCategory.VALID_TOO_LONG.getTag());
            } else {
                setCategoryTag(ReadPairCategory.VALID_PAIR.getTag()); // only if hybrid fragment has the right size it is categorized as valid
            }
        }
    }


    /**
     * @return Linear genomic distance between 5' end positions of the reads.
     */
    private int getDistanceBetweenFivePrimeEnds() {
        if(isTrans()){
            logger.error("Distance between 5' ends on different chromosomes is not defined!");
            return -1;
        } else {
            if(this.getFivePrimeEndPosOfR1() < getFivePrimeEndPosOfR2()) {
                return getFivePrimeEndPosOfR2()-getFivePrimeEndPosOfR1();
            } else {
                return getFivePrimeEndPosOfR1()-getFivePrimeEndPosOfR2();
            }
        }
    }

    public Integer getForwardDigestStart() {return this.digestPair.forward().getStartpos();}
    public Integer getForwardDigestEnd() {return this.digestPair.forward().getEndpos();}
    public boolean forwardDigestIsActive() {return this.digestPair.forward().isActive();}

    public Integer getReverseDigestStart() {return this.digestPair.reverse().getStartpos();}
    public Integer getReverseDigestEnd() {return this.digestPair.reverse().getEndpos();}
    public boolean reverseDigestIsActive() {return this.digestPair.reverse().isActive();}
}
