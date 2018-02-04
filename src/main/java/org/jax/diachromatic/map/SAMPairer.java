package org.jax.diachromatic.map;


import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileWriterFactory;


import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Log;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;
import org.jax.diachromatic.io.Commandline;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * This class takes as input two SAM files that have been created by {@code bowtie2} from the
 * truncated FASTQ files produced by {@link org.jax.diachromatic.command.TruncateCommand}. Its purpose
 * is to rejoin the pairs of reads that correspond to the chimeric fragments in the input files and
 * to perform Q/C and filtering on the reads to remove characteristic Hi-C artefacts.
 * Note that we have made several of the functions in this class package access for testing purposes
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.2 (2018-01-06)
 */
public class SAMPairer {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(SAMPairer.class);
    /**
     * Version of diachromatic. This is initialized within the command line class on the basis of the program
     * version given in the pom.xml file. A default number of zero is given in case initialization doesnt work.
     */
    private static String VERSION = "0.0";
    /**
     * Path to the SAMfile representing the forward read of a paired end experiment. The SAM files should have been
     * processed with the truncate command of this package
     */
    private String samPath1;
    /**
     * Path to the SAMfile representing the reverse read of a paired end experiment. The SAM files should have been
     * processed with the truncate command of this package
     */
    private String samPath2;
    /**
     * Number of unmapped forward reads (SAM flag==4)
     */
    private int n_unmapped_read1 = 0;
    /**
     * Number of unmapped reverse reads (Note: SAM flag is also 4, because the reverse reads are mapped as single-end reads).
     */
    private int n_unmapped_read2 = 0;
    /**
     * Number of pairs with 1 or 2 unmapped reads (these pairs are discarded from further analysis).
     */
    private int n_unmapped_pair = 0;
    /**
     * Number of forward reads that were multimapped (had an XS tag)
     */
    private int n_multimapped_read1 = 0;
    /**
     * Number of reverse reads that were multimapped (had an XS tag)
     */
    private int n_multimapped_read2 = 0;
    /**
     * Number of  read pairs where one or two reads were multimapped (had an XS tag)
     */
    private int n_multimappedPair = 0;

    private int n_could_not_assign_to_digest = 0;
    /** Number of readpairs whose insert was found to have a size above the threshold defined  in {@link ReadPair}.*/
    private int n_insert_too_long = 0;
    /** Number of readpairs whose insert was found to have a size below the threshold defined  in {@link ReadPair}.*/
    private int n_insert_too_short = 0;
    /**
     * Number of circularized reads, a type of artefact where the ends of one fragment ligate with each other.
     */
    private int n_circularized_read = 0;
    /**
     * Number of dangling end reads, a type of artefact where one end of the read is internal and the other is at the
     * end of a restriction fragment.
     */
    private int n_same_dangling_end = 0;

    private int n_same_internal = 0;

    private int n_religation = 0;

    private int n_contiguous = 0;

    private int n_duplicate=0;
    /** Number of reads that pass all quality filters.*/
    private int n_good = 0;
    /**
     * Total number of reads TODO do we mean each read of the paired end reads?
     */
    private int n_total = 0;


    /**
     * Key: chromosome; value: a list of {@link Digest} objects on the chromosome.
     */
    private Map<String, List<Digest>> digestmap = null;
    /**
     * A reader for the forward reads.
     */
    final private SamReader reader1;
    /**
     * A reader for the reverse reads.
     */
    final private SamReader reader2;
    /**
     * Iterator over reads from {@link #reader1}.
     */
    final private Iterator<SAMRecord> it1;
    /**
     * Iterator over reads from {@link #reader2}.
     */
    final private Iterator<SAMRecord> it2;
    /** This will be used to keep a record of valid ditags in order to throw out duplicates. */
    private Set<DiTag> ditagSet;
    /** Count up the number of errors encountered in our reads. THe key is the type of error, and the value is
     * the count over the entire pair of SAM files.
     */
    private Map<QCCode,Integer> errorCounts;


    /**
     * Handle to write valid reads. Used in {@link #inputSAMfiles()}.
     */
    private SAMFileWriter validReadsWriter;
    /**
     * Handle to write invalid reads. Used in {@link #inputSAMfiles()}.
     */
    private SAMFileWriter rejectedReadsWriter;

    private String validBamFileName = "diachromatic.valid.bam";

    private String rejectedBamFileName = "diachromatic.rejected.bam";
    /**If set to true, rejected readpairs are output to {@link #rejectedBamFileName} .*/
    private final boolean outputRejectedReads;

    /**
     * @param sam1    SAM file for the truncated "forward" reads
     * @param sam2    SAM file for the truncated "reverse" reads
     * @param digests see {@link #digestmap}.
     */
    public SAMPairer(String sam1, String sam2, Map<String, List<Digest>> digests, boolean outputRejected) {
        samPath1 = sam1;
        samPath2 = sam2;
        reader1 = SamReaderFactory.makeDefault().open(new File(samPath1));
        reader2 = SamReaderFactory.makeDefault().open(new File(samPath2));
        it1 = reader1.iterator();
        it2 = reader2.iterator();
        digestmap = digests;
        outputRejectedReads = outputRejected;
        VERSION = Commandline.getVersion();
        initializeErrorMap();
    }

    /**
     * An iterator over pairs of SAMRecords -- similar to "next()" in a standard iterator, but will return a pair
     * of SAMRecord objects. Both files must be equally long. This function will return null of there is any issue with
     * with of the individual iterators.
     *
     * @return A {@link ReadPair}, i.e., a pair of SAMRecord objects representing the forward and the reverse reads.
     */
    ReadPair getNextPair() {
        if (it1.hasNext() && it2.hasNext()) {
            SAMRecord record1 = it1.next();
            SAMRecord record2 = it2.next();
            return new ReadPair(record1, record2);
        } else {
            return null;
        }
    }



    /**
     * Input the pair of truncated SAM files.
     * As a side effect, write invalid reads to {@link #rejectedBamFileName}.
     */
    public void inputSAMfiles() throws IOException {

        SAMFileHeader header = reader1.getFileHeader();
        String programGroupId = "@PG\tID:Diachromatic\tPN:Diachromatic\tVN:" + VERSION;
        SAMProgramRecord programRecord = new SAMProgramRecord(programGroupId);
        header.addProgramRecord(programRecord);

        // init BAM outfile
        boolean presorted = false;
        this.validReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(validBamFileName));
        if(outputRejectedReads) {
            this.rejectedReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(rejectedBamFileName));
        }
        final ProgressLogger pl = new ProgressLogger(log, 1000000);

        ReadPair pair;
        while ((pair = getNextPair())!= null) {
            n_total++;
            try {
                // first check whether both reads were mapped.
                if (pair.readPairUniquelyMapped()) {
                    pair.pairReads();
                } else {
                    updateErrorMap(pair.getErrorCodes());
                    continue; // discard this read and go to the next one
                }
                // If we get here, then both reads were uniquely mappable.
                if (is_valid(pair)) {
                    // set the SAM flags to paired-end
                    if (! DiTag.isDuplicate(pair)) { // check for duplicate reads
                        validReadsWriter.addAlignment(pair.forward());
                        validReadsWriter.addAlignment(pair.reverse());
                    } else {
                        n_duplicate++;
                    }
                    n_good++;
                } else {
                    updateErrorMap(pair.getErrorCodes());
                    if (outputRejectedReads) {
                        rejectedReadsWriter.addAlignment(pair.forward());
                        rejectedReadsWriter.addAlignment(pair.reverse());
                    }
                    // discard this read and go to the next one
                }
            } catch (DiachromaticException e) {
                logger.error(e.getMessage()); // todo refactor
            }
        }
        validReadsWriter.close();
        if(outputRejectedReads) {
            rejectedReadsWriter.close();
        }
    }




    /**
     * Decide if a candidate readpair (pair of SAMRecord objects) is valid according to the rules for capture Hi-C
     * @return true if the read is valid
     */
    boolean is_valid(ReadPair readpair) throws DiachromaticException {
        //-1. check whether we can find restriction digests that match the read pair.
        DigestPair digestPair = getDigestPair(readpair);
        if (digestPair == null) {
            readpair.setInvalidDigest();
            return false;
        }
        //-2. Check that calculated insert size is valid
        if (! readpair.hasValidInsertSize(digestPair)) {
            return false;
        }
        //-3. Check if both reads are on the same fragment
        // There are three subclasses of this, all are not valid.
//        if (digestPair.forward().equals(digestPair.reverse())) { // both reads in same restriction fragment.
//            if (readpair.selfLigation()) {
//                return false;
//            } else if (readpair.danglingEnd(digestPair)) {
//                return false;
//            } else {
//                // if we get here, we have reads from the same digest that are not circularized and are not dangling end, so
//                // they must be same_internal
//                readpair.setSameInternal();
//                return false;
//            }
//        }
        if (readpair.bothReadsLocatedOnSameRestrictionFragment(digestPair)) {
            return false;
        }
        // If we get here, then the reads do not map to the same restriction fragment.
        //-4. If the reads map to neighboring fragments, then there may be a religation.
        if (readpair.religation(digestPair)) {
            return false;
        }
        //-5. If we get here, we are on different fragments and the two fragments are not direct neighbors. If they are located
        // within one expected fragment size, then they are contiguous sequences that were not properly digested
        if (readpair.contiguous()) {
            return false;
        }

        // when we get here, we have ruled out artefacts
        return true;
    }











    /**
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

    /** The map {@link #errorCounts} is initialize by setting the counts for all elements to zero. */
    private void initializeErrorMap() {
        this.errorCounts=new HashMap<>();
        for (QCCode ec : QCCode.values()) {
            errorCounts.put(ec,0);
        }
    }

    /**
     * Increment the error code for the errors encounted in a read pair.
     * @param errors Set of errors encountered for some read pair.
     */
    private void updateErrorMap(Set<QCCode> errors) {
        for (QCCode ec : errors) {
            errorCounts.put(ec,1+errorCounts.get(ec));
        }
    }


    public void printStatistics() {
        logger.trace(String.format("n_total pairs=%d\n", n_total));
        logger.trace(String.format("n_unmapped_read1=%d", n_unmapped_read1));
        logger.trace(String.format("n_unmapped_read2=%d", n_unmapped_read2));
        logger.trace(String.format("n_unmapped_pair=%d (%.1f%%)", n_unmapped_pair, (100.0 * n_unmapped_pair / n_total)));

        logger.trace(String.format("n_multimapped_read1=%d", n_multimapped_read1));
        logger.trace(String.format("n_multimapped_read2=%d", n_multimapped_read2));
        logger.trace(String.format("n_multimappedPair=%d (%.1f%%)", n_multimappedPair, (100.0 * n_multimappedPair / n_total)));
        logger.trace(String.format("n_could_not_assign_to_digest=%d (%.1f%%)", n_could_not_assign_to_digest, (100.0 * n_could_not_assign_to_digest / n_total)));
        logger.trace(String.format("n_insert_too_long=%d  (%.1f%%)", n_insert_too_long, (100.0 * n_insert_too_long / n_total)));
        logger.trace(String.format("n_insert_too_short=%d  (%.1f%%)", n_insert_too_short, (100.0 * n_insert_too_short / n_total)));
        logger.trace(String.format("n_circularized_read=%d", n_circularized_read));
        logger.trace(String.format("n_same_dangling_end=%d", n_same_dangling_end));
        logger.trace(String.format("n_same_internal=%d", n_same_internal));
        logger.trace(String.format("n_religation=%d", n_religation));
        logger.trace(String.format("n_contiguous=%d", n_contiguous));
        logger.trace(String.format("n_good=%d (%.1f%%)", n_good, (100.0 * n_good / n_total)));


    }


}
