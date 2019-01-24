package org.jax.diachromatic.align;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.io.Commandline;
import java.io.*;
import java.util.*;

/**
 * This class takes as input two SAM files that have been created by {@code bowtie2} from the
 * truncated FASTQ files produced by {@link org.jax.diachromatic.command.TruncateCommand}. Its purpose
 * is to rejoin the pairs of reads that correspond to the chimeric fragments in the input files and
 * to perform Q/C and filtering on the reads to remove characteristic Hi-C artefacts.
 * Note that we have made several of the functions in this class package access for testing purposes.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.2 (2018-01-06)
 */
public class Aligner {

    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * Version of Diachromatic. This is initialized within the command line class on the basis of the program
     * version given in the pom.xml file. A default number of zero is given in case initialization doesn't work.
     */
    private static String VERSION = "0.0";

    /**
     * Path to the SAM file representing the forward (R1) and reverse (R2) read of a paired end experiment.
     * The SAM files should have been processed with the truncate command of this package.
     */
    private String sam_path_R1;
    private String sam_path_R2;

    /**
     * Number of unmapped forward reads (SAM flag==4 for unmapped). Note: SAM flag is also 4 for the reverse read,
     * because the reads are mapped independently as single-end reads.
     */
    private int n_unmapped_R1 = 0;
    private int n_unmapped_R2 = 0;

    /**
     * Number of read pairs for which at least one read is unmapped.
     */
    private int  n_unmappedPair = 0;

    /**
     * Number of forward reads that were multi-mapped (had an XS tag).
     */
    private int n_multimapped_R1 = 0;
    private int n_multimapped_R2 = 0;

    /**
     * Number of read pairs for which at least one read is multi-mapped.
     */
    private int n_multimappedPair = 0;

    /**
     * Number of paired pairs, i.e. read pairs for which both reads can be mapped uniquely.
     */
    private int n_paired = 0;

    /**
     * Number of unique paired pairs, i.e. same as 'n_paired' but after removal of duplicates.
     */
    private int n_paired_unique = 0;

    /**
     * Number of duplicated read pairs that were removed.
     */
    private int n_paired_duplicated = 0;

    /**
     * Number trans read pairs, i.e. the two reads map to different chromosomes. Trans read pairs are counted after
     * removal of duplicates only.
     */
    private int n_paired_unique_trans = 0;

    /**
     * Total number of truncated read pairs that passed to Diachromatic with the subcommand align.
     */
    private int n_total_input_read_pairs = 0;

    /**
     * Lower and upper bounds for sizes of hybrid fragments. Are passed as arguments to the constructor of {@link ReadPair}
     * and compared to the calculated insert size in order to categorize given read pairs as 'too short' or 'too long'.
     */
    private Integer lowerFragSize;
    private Integer upperFragSize;

    /**
     * Largest calculated insert size represented in the distribution of fragment sizes.
     */
    private static int FRAG_SIZE_LIMIT = 10000;

    /**
     * Arrays that represent size distributions.
     */
    private int[] fragSizesAllPairs =  new int[FRAG_SIZE_LIMIT+1];
    private int[] fragSizesHybridActivePairs =  new int[FRAG_SIZE_LIMIT+1];

    /**
     * If true, multi-mapped reads are defined as those for which a score of a second best hit is reported by bowtie2,
     * i.e. the corresponding SAM record has an XS attribute.
     *
     * If false, a less stringent criterion for uniqueness is used. Reads for which the XS attribute is reported are
     * still categorized as unique, if they have a MAPQ of at least 30 and the difference between AS and XS at least 10.
     */
    private boolean useStringentUniqueSettings;

    /**
     * If true, rejected read pairs are output to an extra BAM file {@link #outputBAMrejected}.
     */
    private final boolean outputRejectedReads;

    /**
     * Filenames (including path) for output BAM files and for text file containing statistics about the alignment and
     * filtering step.
     */
    private String outputBAMvalid, outputBAMrejected, outputTxtStats, outputFragSizesCountsRscript;
    private String filenamePrefix;


    /**
     * Count variables for different read pair categories.
     */
    private int n_un_ligated_pair = 0;
    private int n_self_ligated_pair = 0;

    private int n_valid_pair = 0;
    private int n_valid_too_short = 0;
    private int n_valid_too_long = 0;

    private int n_total_dangling = 0;

    private int n_un_ligated_dangling = 0;
    private int n_self_ligated_dangling = 0;
    private int n_valid_dangling = 0;


    /**
     * Central customized auxiliary class of Diachromatic. Contains information about all restriction fragments of the
     * genome. Is contructed from the digest file produced using GOPHER.
     */
    private  DigestMap digestMap;

    /**
     * HTS-JDK SAM reader objects for R1 and R2.
     */
    final private SamReader sam_reader_R1;
    final private SamReader sam_reader_R2;

    /**
     * HTS-JDK SAM file handles to write valid and rejected read pairs. Used in function {@link #inputSAMfiles()}.
     */
    private SAMFileWriter validReadsWriter;
    private SAMFileWriter rejectedReadsWriter;

    /**
     * Iterator over reads from R1 {@link #sam_reader_R1} and R2 {@link #sam_reader_R2}.
     */
    final private Iterator<SAMRecord> it1;
    final private Iterator<SAMRecord> it2;

    /**
     * @param sam1    SAM file for the truncated "forward" reads
     * @param sam2    SAM file for the truncated "reverse" reads
     */
    public Aligner(String sam1, String sam2, boolean outputRejected, String outputPathPrefix, DigestMap digestMap, Integer lowerFragSize, Integer upperFragSize, String filenamePrefix, boolean useStringentUniqueSettings) {
        sam_path_R1 = sam1;
        sam_path_R2 = sam2;
        sam_reader_R1 = SamReaderFactory.makeDefault().open(new File(sam_path_R1));
        sam_reader_R2 = SamReaderFactory.makeDefault().open(new File(sam_path_R2));
        it1 = sam_reader_R1.iterator();
        it2 = sam_reader_R2.iterator();
        this.digestMap = digestMap;
        outputRejectedReads = outputRejected;
        this.lowerFragSize=lowerFragSize;
        this.upperFragSize=upperFragSize;
        this.filenamePrefix=filenamePrefix;
        this.useStringentUniqueSettings=useStringentUniqueSettings;

        Arrays.fill(fragSizesAllPairs, 0);
        Arrays.fill(fragSizesHybridActivePairs, 0);

        VERSION = Commandline.getVersion();
        createOutputNames(outputPathPrefix);
    }

    /**
     * An iterator over pairs of SAMRecords -- similar to "next()" in a standard iterator, but will return a pair
     * of SAMRecord objects. Both files must be equally long. This function will return null of there is any issue with
     * with of the individual iterators.
     *
     * @return A {@link ReadPair}, i.e., a pair of SAMRecord objects representing the forward and the reverse reads.
     */
    ReadPair getNextPair() throws DiachromaticException {
        if (it1.hasNext() && it2.hasNext()) {
            SAMRecord record1 = it1.next();
            SAMRecord record2 = it2.next();
            return new ReadPair(record1, record2, digestMap, lowerFragSize, upperFragSize, useStringentUniqueSettings);
        } else {
            return null;
        }
    }

    /**
     * Input the pair of truncated SAM files. We will add the PG groups of both
     * SAM files to the header of the output file, and also add a line about the Diachromatic processing.
     * As a side effect, write invalid reads to .
     */
    public void inputSAMfiles() throws IOException, DiachromaticException {

        SAMFileHeader header = sam_reader_R1.getFileHeader();
        SAMFileHeader header2 = sam_reader_R2.getFileHeader();
        // first add program records from the reverse SAM file
        List<SAMProgramRecord> pgList = header2.getProgramRecords();
        for (SAMProgramRecord spr : pgList) {
            //header.addProgramRecord(spr);
        }
        // now add the new program record from Diachromatic
        String programGroupId = "Diachromatic\tPN:Diachromatic\tVN:" + VERSION;
        SAMProgramRecord programRecord = new SAMProgramRecord(programGroupId);
        header.addProgramRecord(programRecord);

        // init BAM outfile
        boolean presorted = false;
        this.validReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMvalid));
        if(outputRejectedReads) {
            this.rejectedReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMrejected));
        }

        DeDupMap dedup_map = new DeDupMap(true);

        ReadPair pair;


        while ((pair = getNextPair())!= null) {

            n_total_input_read_pairs++;

            // first check whether both reads were mapped
            if(pair.isUnMappedR1()) {
                n_unmapped_R1++;}
            if(pair.isUnMappedR2()) {
                n_unmapped_R2++;}
            if(pair.isUnMappedR1()||pair.isUnMappedR2()) {n_unmappedPair++;}
            if(pair.isMultiMappedR1()) {
                n_multimapped_R1++;}
            if(pair.isMultiMappedR2()) {
                n_multimapped_R2++;}
            if(pair.isMultiMappedR1()||pair.isMultiMappedR2()) {n_multimappedPair++;}

            // count categories of pairs
            if(pair.isPaired()) {

                n_paired++;

                // de-duplication starts with paired pairs
                if(dedup_map.hasSeen(pair)) {
                    n_paired_duplicated++;
                    continue;
                }

                n_paired_unique++;

                if(pair.getCategoryTag().equals("VP")) {
                    n_valid_pair++;}
                if(pair.getCategoryTag().equals("UL")) {n_un_ligated_pair++;}
                if(pair.getCategoryTag().equals("SL")) {n_self_ligated_pair++;}
                if(pair.getCategoryTag().equals("TS")) {
                    n_valid_too_short++;}
                if(pair.getCategoryTag().equals("TL")) {
                    n_valid_too_long++;}

                if(pair.isDanglingEnd()) {
                    n_total_dangling++;
                    if(pair.getCategoryTag().equals("VP")) {
                        n_valid_dangling++;}
                    if(pair.getCategoryTag().equals("UL")) {
                        n_un_ligated_dangling++;}
                    if(pair.getCategoryTag().equals("SL")) {
                        n_self_ligated_dangling++;}
                }

                if(pair.isTrans()) {
                    n_paired_unique_trans++;
                }
            }

            // Both reads were uniquely mapped, otherwise the this pair is not paired. If so, continue.
            if(!pair.isPaired()) {continue;}


            // count sizes of all fragments
            Integer incrementFragSize = pair.getCalculatedInsertSize();
            if(FRAG_SIZE_LIMIT<incrementFragSize) { incrementFragSize = FRAG_SIZE_LIMIT; }
            fragSizesAllPairs[incrementFragSize]++;

            // count sizes of all hybrid active fragments
            if((pair.forwardDigestIsActive() & !pair.reverseDigestIsActive()) || (!pair.forwardDigestIsActive() & pair.reverseDigestIsActive())) {
                fragSizesHybridActivePairs[incrementFragSize]++;
            }

            if(pair.getCategoryTag().equals("VP")){

                validReadsWriter.addAlignment(pair.forward());
                validReadsWriter.addAlignment(pair.reverse());
           } else {
                if (outputRejectedReads) {
                    rejectedReadsWriter.addAlignment(pair.forward());
                    rejectedReadsWriter.addAlignment(pair.reverse());
                }
            }
        }

        validReadsWriter.close();
        if(outputRejectedReads) {
            rejectedReadsWriter.close();
        }

        printFragmentLengthDistributionRscript(fragSizesAllPairs, fragSizesHybridActivePairs);

        //dedup_map.printDeDupStatistics(n_paired_duplicated);
    }

    private void printFragmentLengthDistributionRscript(int[] fragSizesAllPairs, int[] fragSizesHybridActivePairs ) throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputFragSizesCountsRscript));

        printStream.print("length<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT-1; i++) {
            printStream.print(i + ",");
        }
        printStream.print(FRAG_SIZE_LIMIT-1 + ")\n");

        printStream.print("fragSizesAllPairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT-1; i++) {
            printStream.print(fragSizesAllPairs[i] + ",");
        }
        printStream.print(fragSizesAllPairs[FRAG_SIZE_LIMIT-1] + ")\n");

        printStream.print("fragSizesHybridActivePairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT-1; i++) {
            printStream.print(fragSizesHybridActivePairs[i] + ",");
        }
        printStream.print(fragSizesHybridActivePairs[FRAG_SIZE_LIMIT-1] + ")\n");

        printStream.print("\n");
        printStream.print("cairo_pdf(\"");
        printStream.print(filenamePrefix);
        printStream.print(".pdf\")\n");

        printStream.print("MAIN=\"");
        printStream.print(filenamePrefix);
        printStream.print("\"\n");

        printStream.print("XLIM<-c(0,1000)\n");

        printStream.print("YLIM<-max(max(fragSizesAllPairs[10:1000]),max(fragSizesHybridActivePairs[10:1000]))\n");

        printStream.print("plot(length, fragSizesAllPairs, xlim=XLIM, type=\"l\", ylim=c(0,YLIM), ylab=NA, xlab=NA, axes=FALSE)\n");

        printStream.print("par(new=TRUE)\n");

        printStream.print("plot(length,fragSizesHybridActivePairs,main=MAIN, xlim=XLIM,type=\"l\", ylim=c(0,YLIM),col=\"red\",xlab=\"Size (nt)\",ylab=\"Fragment count\")\n");

        printStream.print("PREDOM_FRAG_SIZE<-which(max(fragSizesAllPairs)==fragSizesAllPairs)\n");
        printStream.print("abline(v=PREDOM_FRAG_SIZE)\n");

        printStream.print("PREDOM_ACTIVE_FRAG_SIZE<-numeric()\n");
        printStream.print("if(0<sum(fragSizesHybridActivePairs)) {\n");
        printStream.print("PREDOM_ACTIVE_FRAG_SIZE<-which(max(fragSizesHybridActivePairs)==fragSizesHybridActivePairs)\n");
        printStream.print("abline(v=PREDOM_ACTIVE_FRAG_SIZE)\n");
        printStream.print("} else {PREDOM_ACTIVE_FRAG_SIZE <-0}\n");

        printStream.print("LEGEND_ALL<-paste(\"All fragments (\",PREDOM_FRAG_SIZE,\")\",sep=\"\")\n");
        printStream.print("LEGEND_HYBRID_ACTIVE<-paste(\"Hybrid active fragments (\",PREDOM_ACTIVE_FRAG_SIZE,\")\",sep=\"\")\n");
        printStream.print("legend(\"topright\",legend=c(LEGEND_ALL, LEGEND_HYBRID_ACTIVE), col=c(\"black\", \"red\"), lty=1, bg = \"white\")\n");

        printStream.print("dev.off()\n");
    }


    public void printStatistics() throws FileNotFoundException {

        logger.trace(String.format("n_total_input_read_pairs pairs=%d\n", n_total_input_read_pairs));

        logger.trace(String.format("n_unmapped_R1=%d", n_unmapped_R1));
        logger.trace(String.format("n_unmapped_R2=%d\n", n_unmapped_R2));

        logger.trace(String.format("n_multimapped_R1=%d", n_multimapped_R1));
        logger.trace(String.format("n_multimapped_R2=%d\n", n_multimapped_R2));
        logger.trace(String.format("n_multimappedPair=%d\n", n_multimappedPair));

        logger.trace(String.format("n_paired=%d (%.1f%%)\n", n_paired, (100.0 * n_paired / n_total_input_read_pairs)));

        logger.trace(String.format("n_total_dangling=%d", n_total_dangling));



        logger.trace(String.format("n_valid_too_long=%d  (%.1f%%)", n_valid_too_long, (100.0 * n_valid_too_long / n_paired_unique)));
        logger.trace(String.format("n_valid_too_short=%d  (%.1f%%)\n", n_valid_too_short, (100.0 * n_valid_too_short / n_paired_unique)));

        logger.trace("----\n");
        logger.trace(String.format("n_valid_pair=%d", n_valid_pair));
        logger.trace(String.format("n_un_ligated_pair=%d", n_un_ligated_pair));
        logger.trace(String.format("n_self_ligated_pair=%d", n_self_ligated_pair));
       logger.trace(String.format("n_valid_too_short=%d", n_valid_too_short));
        logger.trace(String.format("n_valid_too_long=%d", n_valid_too_long));
        logger.trace(String.format("n_total_input_read_pairs=%d", n_valid_pair + n_valid_too_short + n_valid_too_long +n_un_ligated_pair+n_self_ligated_pair));
        logger.trace("----\n");
        logger.trace(String.format("n_total_dangling=%d", n_total_dangling));
        logger.trace(String.format("n_valid_dangling=%d", n_valid_dangling));
        logger.trace(String.format("n_un_ligated_dangling=%d", n_un_ligated_dangling));
        logger.trace(String.format("n_self_ligated_dangling=%d", n_self_ligated_dangling));

        logger.trace("----\n");

        logger.trace(String.format("n_valid_pair=%d (%.1f%%)", n_valid_pair, (100.0 * n_valid_pair /n_paired_unique)));
        logger.trace("");
        logger.trace("Total number of pairs: " + (n_un_ligated_pair+n_self_ligated_pair+n_valid_too_short+n_valid_too_long+n_valid_pair));
        logger.trace("");
        logger.trace("\t" + "Enrichment Coefficients:");
        logger.trace("\t\t" + "Yield of Valid Pairs (YVP): " + String.format("%.2f%%", 100.0* n_valid_pair / n_total_input_read_pairs));
        logger.trace("\t\t" + "Cross-ligation coefficient (CLC): " + String.format("%.2f%%", 100.0* n_paired_unique_trans /n_paired_unique));
        logger.trace("\t\t" + "Re-ligation coefficient (RLC): " + String.format("%.2f%%", 100.0*(n_paired_unique- n_total_dangling)/n_paired_unique));
        logger.trace("\t\t" + "Pair duplication rate: " + String.format("%.2f%%", 100.0*n_paired_duplicated/n_paired));
        logger.trace("n_duplicate: " + n_paired_duplicated);
        logger.trace("");

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTxtStats));

        printStream.print("Summary statistics\n");
        printStream.print("==================\n\n");
        printStream.print("\n");
        printStream.print("Alignment statistics\n");
        printStream.print("--------------------\n");
        printStream.print("\n");
        printStream.print("Total number of read pairs processed:\t" + n_total_input_read_pairs + "\n");

        printStream.print("Number of unmapped read pairs:\t" + n_unmappedPair + String.format(" (%.2f%%)", 100.0*n_unmappedPair/ n_total_input_read_pairs) + "\n");
        printStream.print("\tNumber of unpapped R1 reads:\t" + n_unmapped_R1 + "\n");
        printStream.print("\tNumber of unpapped R2 reads:\t" + n_unmapped_R2 + "\n");

        printStream.print("Number of multimapped read pairs:\t" + n_multimappedPair + String.format(" (%.2f%%)", 100.0*n_multimappedPair/ n_total_input_read_pairs) + "\n");
        printStream.print("\tNumber of multimapped R1 reads:\t" + n_multimapped_R1 + "\n");
        printStream.print("\tNumber of multimapped R2 reads:\t" + n_multimapped_R2 + "\n");
        printStream.print("Note:\tThere may an overlap between unmapped and multimapped pairs." + "\n");

        printStream.print("Number of paired read pairs:\t" + n_paired + String.format(" (%.2f%%)", 100.0*n_paired/ n_total_input_read_pairs) + "\n");
        printStream.print("\tNumber of unique paired read pairs:\t" + n_paired_unique + "\n");
        printStream.print("\tNumber of duplicated pairs:\t" + n_paired_duplicated + "\n");

        printStream.print("\n");
        printStream.print("Artifact statistics\n");
        printStream.print("-------------------\n");
        printStream.print("\n");

        printStream.print("Dangling end total:\t" + n_total_dangling + String.format(" (%.2f%%)", 100.0* n_total_dangling /n_paired_unique) + "\n");
        printStream.print("Note: Dangling ends may occur in all categories, and a read pair with a dangling end can still be valid." + "\n\n");


        printStream.print("Un-ligated:\t" + n_un_ligated_pair + String.format(" (%.2f%%)", 100.0*n_un_ligated_pair/n_paired_unique) + "\n");
        printStream.print("Self-ligated:\t" + n_self_ligated_pair + String.format(" (%.2f%%)", 100.0*n_self_ligated_pair/n_paired_unique) + "\n");
        printStream.print("Calculated size of hybrid fragment too short:\t" + n_valid_too_short + String.format(" (%.2f%%)", 100.0* n_valid_too_short /n_paired_unique) + "\n");
        printStream.print("Calculated size of hybrid fragment too long:\t" + n_valid_too_long + String.format(" (%.2f%%)", 100.0* n_valid_too_long /n_paired_unique) + "\n");
        printStream.print("Valid:\t" + n_valid_pair + String.format(" (%.2f%%)", 100.0* n_valid_pair /n_paired_unique) + "\n");
        printStream.print("Note: These five categories are disjoint subsets of the unique paired read pairs." + "\n\n");

        printStream.print("\n");
        printStream.print("Quality metrics for experimental trouble shooting \n");
        printStream.print("--------------------------------------------------\n");
        printStream.print("\n");
        printStream.print("Yield of Valid Pairs (YVP):\t" + String.format("%.2f%%", 100.0* n_valid_pair / n_total_input_read_pairs) + "\n");
        printStream.print("Cross-ligation coefficient (CLC):\t" + String.format("%.2f%%", 100.0* n_paired_unique_trans /n_paired_unique) + "\n");
        printStream.print("Re-ligation coefficient (RLC):\t" + String.format("%.2f%%", 100.0*(n_paired_unique-n_total_dangling)/n_paired_unique) + "\n");
        printStream.print("Hi-C pair Duplication Rate (HPDR):\t" + String.format("%.2f%%", 100.0*n_paired_duplicated/n_paired) + "\n");
        printStream.print("\n");
    }

    private void createOutputNames(String outputPathPrefix) {
        outputBAMvalid = String.format("%s.%s", outputPathPrefix, "valid_pairs.aligned.bam");
        outputBAMrejected = String.format("%s.%s", outputPathPrefix, "rejected_pairs.aligned.bam");
        outputFragSizesCountsRscript = String.format("%s.%s", outputPathPrefix, "frag.sizes.counts.script.R");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "align.stats.txt");
    }
}
