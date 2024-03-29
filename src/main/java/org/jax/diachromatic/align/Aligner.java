package org.jax.diachromatic.align;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;

import org.jax.diachromatic.Diachromatic;
import org.jax.diachromatic.exception.DiachromaticException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

/**
 * This class takes as input two SAM files that have been created by {@code bowtie2} from the truncated FASTQ files
 * produced by {@link org.jax.diachromatic.command.TruncateCommand}. Its purpose is to rejoin the pairs of reads that
 * correspond to the chimeric fragments in the input files and to perform Q/C and filtering on the reads to remove
 * characteristic Hi-C artifacts.
 *
 * Note that we have made several of the functions in this class package access for testing purposes.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.2 (2018-01-06)
 */
public class Aligner {
    private static final Logger logger = LoggerFactory.getLogger(Aligner.class);
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * Version of Diachromatic. This is initialized within the command line class on the basis of the program
     * version given in the pom.xml file. A default number of zero is given in case initialization doesn't work.
     */
    private static String VERSION = "0.0";

    /**
     * Total number of truncated read pairs that passed to Diachromatic with the subcommand align.
     */
    private int n_total_input_read_pairs = 0;

    /**
     * Numbers of unmapped forward and reverse reads (SAM flag==4 for unmapped). Note: SAM flag is also 4 for the
     * reverse read, because the reads are mapped independently as single-end reads.
     */
    private int n_unmapped_R1 = 0;
    private int n_unmapped_R2 = 0;

    /**
     * Number of read pairs for which at least one read is unmapped.
     */
    private int  n_unmappedPair = 0;

    /**
     * Numbers of forward and reverse reads that were multi-mapped (had an XS tag).
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
     * Count variables for disjoint read pair categories (see documentation on read the docs).
     */
    private int n_paired_unique_un_ligated = 0;
    private int n_paired_unique_un_ligated_same_internal = 0;
    private int n_paired_unique_self_ligated = 0;
    private int n_paired_unique_self_ligated_same_internal = 0;
    private int n_paired_unique_too_short = 0;
    private int n_paired_unique_too_long = 0;
    private int n_paired_unique_valid = 0;
    private int n_paired_strange_internal = 0;

    /**
     * Number trans read pairs, i.e. the two reads of a given pair map to different chromosomes. Trans read pairs are
     * counted after removal of duplicates only. Trans read pairs cannot be un-ligated or self-ligated by definition.
     */
    private int n_paired_unique_trans = 0;

    /**
     * Number of dangling end read pairs. A read pairs is categorized as dangling end pair if the 5' end position of at
     * least one of the two reads occurs at a distance of at most DANGLING_THRESHOLD = 7 from the next restriction
     * enzyme cutting site. Dangling end read pairs are counted after removal of duplicates only. Dangling end read
     * pairs may occur in all read pair categories.
     */
    private int n_paired_unique_dangling = 0;

    /**
     * Additional experimental count variable for more detailed characterization and sanity checks.
     */
    private int n_paired_unique_un_ligated_dangling = 0;
    private int n_paired_unique_self_ligated_dangling = 0;
    private int n_paired_unique_too_short_dangling = 0;
    private int n_paired_unique_too_long_dangling = 0;
    private int n_paired_unique_valid_dangling = 0;
    private int n_paired_strange_internal_dangling = 0;

    private int n_paired_unique_un_ligated_trans = 0; // should never be incremented
    private int n_paired_unique_self_ligated_trans = 0; // should never be incremented
    private int n_paired_unique_too_short_trans = 0;
    private int n_paired_unique_too_long_trans = 0;
    private int n_paired_unique_valid_trans = 0;
    private int n_paired_strange_internal_trans = 0;

    /**
     * Lower and upper bounds for sizes of chimeric fragments. Are passed as arguments to the constructor of {@link ReadPair}
     * and compared to the calculated insert size in order to categorize given read pairs as 'too short' or 'too long'.
     */
    private Integer lowerFragSize;
    private Integer upperFragSize;

    /**
     * If true, multi-mapped reads are defined as those for which a score of a second best hit is reported by bowtie2,
     * i.e. the corresponding SAM record has an XS attribute.
     *
     * If false, a less stringent criterion for uniqueness is used. Reads for which the XS attribute is reported are
     * still categorized as unique, if they have a MAPQ of at least 30 and the difference between AS and XS at least 10.
     */
    private boolean useStringentUniqueSettings;

    private boolean useRelativeOrientationForDuplicateRemoval = false;

    /**
     * If true, rejected read pairs are summarize to an extra BAM file {@link #outputBAMrejected}.
     */
    private final boolean outputRejectedReads;

    /**
     * Filenames (including path) for summarize BAM files and for text file containing statistics about the alignment and
     * filtering step.
     */
    private String outputBAMvalid, outputBAMrejected, outputTxtStats, outputFragSizesCountsRscript;
    private String filenamePrefix;

    /**
     * Central customized auxiliary class of Diachromatic. Contains information about all restriction fragments of the
     * genome. Is constructed from the digest file produced using GOPHER.
     */
    private  DigestMap digestMap;

    /**
     * Largest calculated insert size represented in the distribution of fragment sizes.
     */
    private static int FRAG_SIZE_LIMIT = 30000;

    /**
     * Arrays that represent size distributions.
     */
    private int[] fragSizesChimericPairs =  new int[FRAG_SIZE_LIMIT+1];
    private int[] fragSizesActiveChimericPairs =  new int[FRAG_SIZE_LIMIT+1];
    private int[] fragSizesUnLigatedPairs =  new int[FRAG_SIZE_LIMIT+1];
    private int[] fragSizesSelfLigatedSameInternalPairs =  new int[FRAG_SIZE_LIMIT+1];

    /**
     * HasMap for Trans/Cis ratio
     */
    Map<String, Integer> cisCounts;
    Map<String, Integer> transCounts;


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
     * Path to the SAM file representing the forward (R1) and reverse (R2) read of a paired end experiment.
     * The SAM files should have been processed with the truncate command of this package.
     */
    private String sam_path_R1;
    private String sam_path_R2;

    /**
     * Constructor of this class.
     *
     * @param sam1 SAM file for the truncated R1 reads
     * @param sam2 SAM file for the truncated R2 reads
     * @param outputRejected If true, an additional BAM file for rejected reads will be created.
     * @param outputPathPrefix Path for summarize including path and file prefix.
     * @param digestMap Custom class of Diachromatic containing information about all digests of the genome.
     * @param lowerFragSize Lower threshold for fragments sizes consistent with sonication parameters.
     * @param upperFragSize Upper threshold for fragments sizes consistent with sonication parameters.
     * @param filenamePrefix Prefix for names of created files.
     * @param useStringentUniqueSettings Use the more stringent definition of multi-mapped reads.
     */
    public Aligner(String sam1, String sam2, boolean outputRejected, String outputPathPrefix, DigestMap digestMap,
                   Integer lowerFragSize, Integer upperFragSize, Integer upperSelfLigationSize, String filenamePrefix,
                   boolean useStringentUniqueSettings) {
        this.sam_path_R1 = sam1;
        this.sam_path_R2 = sam2;
        this.sam_reader_R1 = SamReaderFactory.makeDefault().open(new File(sam_path_R1));
        this.sam_reader_R2 = SamReaderFactory.makeDefault().open(new File(sam_path_R2));
        this.it1 = sam_reader_R1.iterator();
        this.it2 = sam_reader_R2.iterator();
        this.digestMap = digestMap;
        this.outputRejectedReads = outputRejected;
        this.lowerFragSize = lowerFragSize;
        this.upperFragSize = upperFragSize;
        this.filenamePrefix = filenamePrefix;
        this.useStringentUniqueSettings = useStringentUniqueSettings;
        this.useRelativeOrientationForDuplicateRemoval = false;
        ReadPair.setLengthThresholds(lowerFragSize,upperFragSize,upperSelfLigationSize);
        Arrays.fill(fragSizesChimericPairs, 0);
        Arrays.fill(fragSizesActiveChimericPairs, 0);
        Arrays.fill(fragSizesUnLigatedPairs, 0);
        Arrays.fill(fragSizesSelfLigatedSameInternalPairs, 0);



        VERSION = Diachromatic.getVersion();
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
            return new ReadPair(record1, record2, digestMap, useStringentUniqueSettings);
        } else {
            return null;
        }
    }

    /**
     * Input the pair of truncated SAM files. We will add the PG groups of both
     * SAM files to the header of the summarize file, and also add a line about the Diachromatic processing.
     *
     * @throws IOException Required because of generation R script for fragment size distribution.
     * @throws DiachromaticException Required because class ReadPair is used.
     */
    public void inputSAMfiles() throws IOException, DiachromaticException {

        SAMFileHeader header = sam_reader_R1.getFileHeader();
        SAMFileHeader header2 = sam_reader_R2.getFileHeader();

        // add the new program record from Diachromatic
        String programGroupId = "Diachromatic\tPN:Diachromatic\tVN:" + VERSION;
        SAMProgramRecord programRecord = new SAMProgramRecord(programGroupId);
        header.addProgramRecord(programRecord);

        // init BAM outfile
        boolean presorted = false;
        this.validReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMvalid));
        if(outputRejectedReads) {
            this.rejectedReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMrejected));
        }

        DeDupMap dedup_map = new DeDupMap(useRelativeOrientationForDuplicateRemoval);
        ReadPair pair;
        cisCounts = new HashMap<>();
        transCounts = new HashMap<>();

        while ((pair = getNextPair())!= null) {

            n_total_input_read_pairs++;

            if(n_total_input_read_pairs%1000000==0) {
                logger.trace("n_total_input_read_pairs: " + n_total_input_read_pairs);
            }

            if(dedup_map.getNumOfInsertions()%1000000==0 && 0<dedup_map.getNumOfInsertions()) {
                logger.trace("dedup_map.getNumOfInsertions(): " + dedup_map.getNumOfInsertions());
            }

            // first check whether both reads were mapped uniquely
            if (pair.isUnMappedR1()) {
                n_unmapped_R1++;
            }
            if (pair.isUnMappedR2()) {
                n_unmapped_R2++;
            }
            if (pair.isUnMappedR1() || pair.isUnMappedR2()) {
                n_unmappedPair++;
            }
            if (pair.isMultiMappedR1()) {
                n_multimapped_R1++;
            }
            if (pair.isMultiMappedR2()) {
                n_multimapped_R2++;
            }
            if (pair.isMultiMappedR1() || pair.isMultiMappedR2()) {
                n_multimappedPair++;
            }

            // Note: Read pairs with unmapped or multi-mapped reads remain unpaired

            // count categories of paired pairs
            if(pair.isPaired()) {

                n_paired++;

                // de-duplication starts with paired pairs
                if(dedup_map.hasSeen(pair)) {
                    n_paired_duplicated++;
                    continue;
                }

                n_paired_unique++;

                if(pair.getCategoryTag().equals("VP")) {
                    n_paired_unique_valid++;}
                if(pair.getCategoryTag().equals("UL")) {
                    n_paired_unique_un_ligated++;}
                if(pair.getCategoryTag().equals("ULSI")) {
                    n_paired_unique_un_ligated_same_internal++;}
                if(pair.getCategoryTag().equals("SL")) {
                    n_paired_unique_self_ligated++;}
                if(pair.getCategoryTag().equals("SLSI")) {
                    n_paired_unique_self_ligated_same_internal++;}
                if(pair.getCategoryTag().equals("TS")) {
                    n_paired_unique_too_short++;}
                if(pair.getCategoryTag().equals("TL")) {
                    n_paired_unique_too_long++;}
                if(pair.getCategoryTag().equals("SI")) {
                    n_paired_strange_internal++;}

                if(pair.isDanglingEnd()) {
                    n_paired_unique_dangling++;
                    if(pair.getCategoryTag().equals("VP")) {
                        n_paired_unique_valid_dangling++;}
                    if(pair.getCategoryTag().equals("UL")) {
                        n_paired_unique_un_ligated_dangling++;}
                    if(pair.getCategoryTag().equals("SL")) {
                        n_paired_unique_self_ligated_dangling++;}
                    if(pair.getCategoryTag().equals("TS")) {
                        n_paired_unique_too_short_dangling++;}
                    if(pair.getCategoryTag().equals("TL")) {
                        n_paired_unique_too_long_dangling++;}
                    if(pair.getCategoryTag().equals("SI")) {
                        n_paired_strange_internal_dangling++;}
                }

                if(pair.isTrans()) {
                    n_paired_unique_trans++;
                    if(pair.getCategoryTag().equals("VP")) {
                        n_paired_unique_valid_trans++;}
                    if(pair.getCategoryTag().equals("UL")) {
                        n_paired_unique_un_ligated_trans++;}
                    if(pair.getCategoryTag().equals("SL")) {
                        n_paired_unique_self_ligated_trans++;}
                    if(pair.getCategoryTag().equals("TS")) {
                        n_paired_unique_too_short_trans++;}
                    if(pair.getCategoryTag().equals("TL")) {
                        n_paired_unique_too_long_trans++;}
                    if(pair.getCategoryTag().equals("SI")) {
                        n_paired_strange_internal_trans++;}
                    if(pair.getCategoryTag().equals("VP")) { // count trans/cis for chromosome-wise CLC
                        if(transCounts.containsKey(pair.getReferenceSequenceOfR1())) {
                            transCounts.put(pair.getReferenceSequenceOfR1(),transCounts.get(pair.getReferenceSequenceOfR1())+1);
                        } else {
                            transCounts.put(pair.getReferenceSequenceOfR1(),1);
                            if(!cisCounts.containsKey(pair.getReferenceSequenceOfR1())) {
                                cisCounts.put(pair.getReferenceSequenceOfR1(),0);
                            }
                        }
                        if(transCounts.containsKey(pair.getReferenceSequenceOfR2())) {
                            transCounts.put(pair.getReferenceSequenceOfR2(),transCounts.get(pair.getReferenceSequenceOfR2())+1);
                        } else {
                            transCounts.put(pair.getReferenceSequenceOfR2(),1);
                            if(!cisCounts.containsKey(pair.getReferenceSequenceOfR2())) {
                                cisCounts.put(pair.getReferenceSequenceOfR2(),0);
                            }
                        }
                    }
                } else {
                    if(pair.getCategoryTag().equals("VP")) {
                        if (cisCounts.containsKey(pair.getReferenceSequenceOfR1())) {
                            cisCounts.put(pair.getReferenceSequenceOfR1(), cisCounts.get(pair.getReferenceSequenceOfR1()) + 2);
                        } else {
                            cisCounts.put(pair.getReferenceSequenceOfR1(), 2);
                            if(!transCounts.containsKey(pair.getReferenceSequenceOfR1())) {
                                transCounts.put(pair.getReferenceSequenceOfR1(),0);
                            }
                        }
                    }
                }
            } else {
                continue;
            }

            // count sizes of all chimeric fragments including valid, too short and too long
            Integer incrementFragSize = pair.getChimericFragmentSize();
            if(FRAG_SIZE_LIMIT<incrementFragSize) { incrementFragSize = FRAG_SIZE_LIMIT; }
            if(pair.getCategoryTag().equals("VP")||pair.getCategoryTag().equals("TS")||pair.getCategoryTag().equals("TL"))   {

                fragSizesChimericPairs[incrementFragSize]++;

                // count sizes of all active chimeric fragments
                if((pair.forwardDigestIsActive() & !pair.reverseDigestIsActive()) || (!pair.forwardDigestIsActive() & pair.reverseDigestIsActive())) {
                    fragSizesActiveChimericPairs[incrementFragSize]++;
                }
            }

            // count sizes of potentially un-ligated fragments (don't use thresholds to avoid circular argument)
            if(pair.isInwardFacing() && !pair.isTrans()){
                incrementFragSize=pair.getDistanceBetweenFivePrimeEnds();
                if(FRAG_SIZE_LIMIT<incrementFragSize) { incrementFragSize = FRAG_SIZE_LIMIT; }
                fragSizesUnLigatedPairs[incrementFragSize]++;
            }

            // count sizes of potentially un-ligated fragments (don't use thresholds to avoid circular argument)
            if(pair.getCategoryTag().equals("SLSI")){
                incrementFragSize=pair.getSelfLigationFragmentSize();
                if(FRAG_SIZE_LIMIT<incrementFragSize) { incrementFragSize = FRAG_SIZE_LIMIT; }
                fragSizesSelfLigatedSameInternalPairs[incrementFragSize]++;
            }


            // write pair to BAM file
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

        printFragmentLengthDistributionRscript(fragSizesChimericPairs, fragSizesActiveChimericPairs, fragSizesUnLigatedPairs, fragSizesSelfLigatedSameInternalPairs);
        //dedup_map.printDeDupStatistics(n_paired_duplicated);
    }

    /**
     * This function generates an R script that can be used to create a pdf of the distribution of fragments sizes.
     * The sizes for read pairs that belong to selected/active fragments are passed and plotted separately.
     * The purpose of this is to investigate the relationship of fragment size and enrichment.
     *
     * @param fragSizesAllPairs integer array representing the distribution of sizes, e.g. fragSizesChimericPairs[181] corrsponds to the number of fragments of size 180
     * @param fragSizesChimericActivePairs same as fragSizesChimericPairs but only for read pairs for which at least one read maps to an selected/active fragment
     * @throws FileNotFoundException
     */
    private void printFragmentLengthDistributionRscript(int[] fragSizesAllPairs, int[] fragSizesChimericActivePairs, int[] fragSizesUnLigatedPairs, int[] fragSizesSelfLigatedSameInternalPairs) throws FileNotFoundException {

        PrintStream printStream = new PrintStream(new FileOutputStream(outputFragSizesCountsRscript));

        printStream.print("length<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            printStream.print(i + ",");
        }
        printStream.print(FRAG_SIZE_LIMIT + ")\n");

        printStream.print("fragSizesChimericPairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            printStream.print(fragSizesAllPairs[i] + ",");
        }
        printStream.print(fragSizesAllPairs[FRAG_SIZE_LIMIT] + ")\n");

        printStream.print("fragSizesActiveChimericPairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            printStream.print(fragSizesChimericActivePairs[i] + ",");
        }
        printStream.print(fragSizesChimericActivePairs[FRAG_SIZE_LIMIT] + ")\n");

        printStream.print("fragSizesUnLigatedPairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            printStream.print(fragSizesUnLigatedPairs[i] + ",");
        }
        printStream.print(fragSizesUnLigatedPairs[FRAG_SIZE_LIMIT] + ")\n");

        printStream.print("fragSizesSelfLigatedSameInternalPairs<-c(");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            printStream.print(fragSizesSelfLigatedSameInternalPairs[i] + ",");
        }
        printStream.print(fragSizesSelfLigatedSameInternalPairs[FRAG_SIZE_LIMIT] + ")\n");


        printStream.print("\n");
        printStream.print("cairo_pdf(\"");
        printStream.print(filenamePrefix);
        printStream.print(".pdf\", height=6, width=12)\n");
        printStream.print("par(mfrow=c(1,2))\n");

        printStream.print("FRAG_SIZE_LIMIT=");
        printStream.print(FRAG_SIZE_LIMIT+1);
        printStream.print("\n");

        printStream.print("MAIN=\"");
        printStream.print(filenamePrefix);
        printStream.print("\"\n");

        printStream.print("XLIM<-c(0,1000)\n");

        printStream.print("YLIM<-max(max(fragSizesChimericPairs[10:1000]),max(fragSizesActiveChimericPairs[10:1000]))\n");

        printStream.print("plot(length, fragSizesChimericPairs, xlim=XLIM, type=\"l\", ylim=c(0,YLIM), ylab=NA, xlab=NA, axes=FALSE)\n");

        printStream.print("par(new=TRUE)\n");

        printStream.print("plot(length, fragSizesUnLigatedPairs, xlim=XLIM, type=\"l\", ylim=c(0,YLIM), ylab=NA, xlab=NA, axes=FALSE, col=\"blue\")\n");

        printStream.print("par(new=TRUE)\n");

        printStream.print("plot(length,fragSizesActiveChimericPairs,main=\"Size distribution of chimeric and un-ligated fragments\", xlim=XLIM,type=\"l\", ylim=c(0,YLIM),col=\"red\",xlab=\"Size (nt)\",ylab=\"Fragment count\")\n");

        printStream.print("PREDOM_FRAG_SIZE<-which(max(fragSizesChimericPairs)==fragSizesChimericPairs)-1\n");
        //printStream.print("abline(v=PREDOM_FRAG_SIZE,col=\"black\")\n");

        printStream.print("PREDOM_UNLIGATED_FRAG_SIZE<-which(max(fragSizesUnLigatedPairs[1:1000])==fragSizesUnLigatedPairs[1:1000])-1\n");
        //printStream.print((printStream.print("abline(v=PREDOM_UNLIGATED_FRAG_SIZE,col=\"blue\")\n");

        printStream.print("PREDOM_ACTIVE_FRAG_SIZE<-numeric()\n");
        printStream.print("if(0<sum(fragSizesActiveChimericPairs)) {\n");
        printStream.print("PREDOM_ACTIVE_FRAG_SIZE<-which(max(fragSizesActiveChimericPairs)==fragSizesActiveChimericPairs)-1\n");
        printStream.print("abline(v=PREDOM_ACTIVE_FRAG_SIZE,col=\"red\")\n");
        printStream.print("} else {PREDOM_ACTIVE_FRAG_SIZE <-0}\n");


        printStream.print("LEGEND_HYBRID<-paste(\"Chimeric fragments (\",PREDOM_FRAG_SIZE,\")\",sep=\"\")\n");
        printStream.print("LEGEND_ACTIVE<-paste(\"Enriched chimeric fragments (\",PREDOM_ACTIVE_FRAG_SIZE,\")\",sep=\"\")\n");
        printStream.print("LEGEND_UNLIGATED<-paste(\"Un-ligated fragments (\",PREDOM_UNLIGATED_FRAG_SIZE,\")\",sep=\"\")\n");

        printStream.print("legend(\"topright\",legend=c(LEGEND_HYBRID, LEGEND_ACTIVE, LEGEND_UNLIGATED), col=c(\"black\", \"red\", \"blue\"), lty=1, bg = \"white\")\n\n");


        printStream.print("plot(length, fragSizesSelfLigatedSameInternalPairs, xlab=\"Size (nt)\", ylab=\"Fragment count\", type=\"l\", xlim=c(0,20000), col=\"black\", main=\"Size distribution of self-ligated same internal fragments\")\n");

        printStream.print("dev.off()\n");
    }

    /**
     * This function prints summary statistics about the alignment step to the file: prefix.align.stats.txt
     *
     * @throws FileNotFoundException required because of FileOutputStream.
     */
    public void printStatistics() throws FileNotFoundException {

        PrintStream printStream = new PrintStream(new FileOutputStream(outputTxtStats));
        
        printStream.print("total_read_pairs_processed:\t" + n_total_input_read_pairs + "\n");

        printStream.print("unmapped_read_pairs:" + n_unmappedPair + String.format(" (%.2f%%)", 100.0*n_unmappedPair/ n_total_input_read_pairs) + "\n");
        printStream.print("unmapped_R1_reads:" + n_unmapped_R1 + "\n");
        printStream.print("unmapped_R2_reads:" + n_unmapped_R2 + "\n");

        printStream.print("multimapped_read_pairs:" + n_multimappedPair + String.format(" (%.2f%%)", 100.0*n_multimappedPair/ n_total_input_read_pairs) + "\n");
        printStream.print("multimapped_R1_reads:" + n_multimapped_R1 + "\n");
        printStream.print("multimapped_R2_reads:" + n_multimapped_R2 + "\n");

        printStream.print("paired_read_pairs:" + n_paired + String.format(" (%.2f%%)", 100.0*n_paired/ n_total_input_read_pairs) + "\n");
        printStream.print("unique_paired_read_pairs:" + n_paired_unique + "\n");
        printStream.print("duplicated_pairs:\t" + n_paired_duplicated + "\n");
        printStream.print("\n");
        printStream.print("Artifact statistics\n");

        int n_paired_unique_un_ligated_total=n_paired_unique_un_ligated+n_paired_unique_un_ligated_same_internal;
        printStream.print("unligated:" + n_paired_unique_un_ligated_total + String.format(" (%.2f%%)", 100.0* n_paired_unique_un_ligated_total /n_paired_unique) + "\n");
        printStream.print("unligated_by_size:" + n_paired_unique_un_ligated + String.format(" (%.2f%%)", 100.0* n_paired_unique_un_ligated /n_paired_unique) + "\n");
        printStream.print("unligated_same_internal:" + n_paired_unique_un_ligated_same_internal + String.format(" (%.2f%%)", 100.0* n_paired_unique_un_ligated_same_internal /n_paired_unique) + "\n");

        int n_paired_unique_self_ligated_total=n_paired_unique_self_ligated+n_paired_unique_self_ligated_same_internal;
        printStream.print("self_ligated:" + n_paired_unique_self_ligated_total + String.format(" (%.2f%%)", 100.0* n_paired_unique_self_ligated_total /n_paired_unique) + "\n");
        printStream.print("self_ligated_by_size:" + n_paired_unique_self_ligated + String.format(" (%.2f%%)", 100.0* n_paired_unique_self_ligated /n_paired_unique) + "\n");
        printStream.print("self_ligated_same_internal:" + n_paired_unique_self_ligated_same_internal + String.format(" (%.2f%%)", 100.0* n_paired_unique_self_ligated_same_internal /n_paired_unique) + "\n");

        int n_chimeric_fragments=n_paired_unique_too_short+n_paired_unique_too_long+n_paired_unique_valid;

        printStream.print("chimeric:" + n_chimeric_fragments + String.format(" (%.2f%%)", 100.0* n_chimeric_fragments /n_paired_unique) + "\n");
        printStream.print("chimeric_short:" + n_paired_unique_too_short + String.format(" (%.2f%%)", 100.0* n_paired_unique_too_short /n_paired_unique) + "\n");
        printStream.print("chimeric_long:" + n_paired_unique_too_long + String.format(" (%.2f%%)", 100.0* n_paired_unique_too_long /n_paired_unique) + "\n");
        printStream.print("chimeric_valid:" + n_paired_unique_valid + String.format(" (%.2f%%)", 100.0* n_paired_unique_valid /n_paired_unique) + "\n");

        printStream.print("strange_internal:" + n_paired_strange_internal + String.format(" (%.2f%%)", 100.0* n_paired_strange_internal /n_paired_unique) + "\n");

        printStream.print("Note: These four categories are disjoint subsets of all unique paired read pairs, and percentages refer to this superset." + "\n\n");

        printStream.print("dangling_end_pairs_total:" + n_paired_unique_dangling + String.format(" (%.2f%%)", 100.0* n_paired_unique_dangling /n_paired_unique) + "\n");
        printStream.print("Note: Dangling end pairs may occur in all categories, and a read pair with a dangling end can still be valid." + "\n\n");

        printStream.print("trans_pairs_total:" + n_paired_unique_trans + String.format(" (%.2f%%)", 100.0* n_paired_unique_trans /n_paired_unique) + "\n");
        printStream.print("Note: Trans pairs cannot occur in the categories un-ligated and self-ligated but all others." + "\n\n");

//        printStream.print("Quality metrics for experimental trouble shooting\n");
//        printStream.print("-------------------------------------------------\n");
//        printStream.print("\n");
        printStream.print("YVP:\t" + String.format("%.2f%%", 100.0* n_paired_unique_valid / n_total_input_read_pairs) + "\n");
        printStream.print("CLC:\t" + String.format("%.2f%%", 100.0* n_paired_unique_valid_trans/n_paired_unique_valid) + "\n");
        double rlc = 1.0 - 1.0*(n_paired_unique_too_short_dangling+n_paired_unique_too_long_dangling+n_paired_unique_valid_dangling)/(n_paired_unique_too_short+n_paired_unique_too_long+n_paired_unique_valid);
        printStream.print("RLC:\t" + String.format("%.2f%%", 100.0*rlc) + "\n");
        printStream.print("HPDR:\t" + String.format("%.2f%%", 100.0*n_paired_duplicated/n_paired) + "\n");
        printStream.print("\n");


        printStream.print("\n");
        printStream.print("More detailed results and sanity checks\n");
        printStream.print("---------------------------------------\n");
        printStream.print("\n");

        printStream.print("Fractions of dangling end pairs:\n");
        printStream.print(String.format("n_paired_unique_un_ligated_dangling:%d (%.2f%% of all unique un-ligated pairs)\n", n_paired_unique_un_ligated_dangling, (100.0 * n_paired_unique_un_ligated_dangling / n_paired_unique_un_ligated)));
        printStream.print(String.format("n_paired_unique_self_ligated_dangling:%d (%.2f%% of all unique self-ligated pairs)\n", n_paired_unique_self_ligated_dangling, (100.0 * n_paired_unique_self_ligated_dangling / n_paired_unique_self_ligated)));
        printStream.print(String.format("n_paired_unique_too_short_dangling:%d (%.2f%% of all unique valid too short pairs)\n", n_paired_unique_too_short_dangling, (100.0 * n_paired_unique_too_short_dangling / n_paired_unique_too_short)));
        printStream.print(String.format("n_paired_unique_too_long_dangling:%d (%.2f%% of all unique valid too long pairs)\n", n_paired_unique_too_long_dangling, (100.0 * n_paired_unique_too_long_dangling / n_paired_unique_too_long)));
        printStream.print(String.format("n_paired_unique_valid_dangling:%d (%.2f%% of all unique valid pairs)\n", n_paired_unique_valid_dangling, (100.0 * n_paired_unique_valid_dangling / n_paired_unique_valid)));
        printStream.print(String.format("n_paired_strange_internal_dangling:%d (%.2f%% of all unique valid pairs)\n", n_paired_strange_internal_dangling, (100.0 * n_paired_strange_internal_dangling / n_paired_unique_valid)));

        printStream.print("\n");
        printStream.print("Fractions of trans pairs:\n");
        printStream.print(String.format("n_paired_unique_un_ligated_trans:%d (%.2f%% of all unique un-ligated pairs)\n", n_paired_unique_un_ligated_trans, (100.0 * n_paired_unique_un_ligated_trans / n_paired_unique_un_ligated)));
        printStream.print(String.format("n_paired_unique_self_ligated_trans:%d (%.2f%% of all unique self-ligated pairs)\n", n_paired_unique_self_ligated_trans, (100.0 * n_paired_unique_self_ligated_trans / n_paired_unique_self_ligated)));
        printStream.print(String.format("n_paired_unique_too_short_trans:%d (%.2f%% of all unique valid too short pairs)\n", n_paired_unique_too_short_trans, (100.0 * n_paired_unique_too_short_trans / n_paired_unique_too_short)));
        printStream.print(String.format("n_paired_unique_too_long_trans:%d (%.2f%% of all unique valid too long pairs)\n", n_paired_unique_too_long_trans, (100.0 * n_paired_unique_too_long_trans / n_paired_unique_too_long)));
        printStream.print(String.format("n_paired_unique_valid_trans:%d (%.2f%% of all unique valid pairs)\n", n_paired_unique_valid_trans, (100.0 * n_paired_unique_valid_trans / n_paired_unique_valid)));
        printStream.print(String.format("n_paired_strange_internal_trans:%d (%.2f%% of all unique valid pairs)\n", n_paired_strange_internal_trans, (100.0 * n_paired_strange_internal_trans / n_paired_unique_valid)));
        printStream.print(String.format("n_total_trans:%d (%.2f%% of all unique paired read pairs)\n", n_paired_unique_trans, (100.0 * n_paired_unique_trans/n_paired_unique)));

        printStream.print("\n");
        printStream.print("chimeric_fragment_size_count_array:");
        for(int i=0; i<1000; i++) {
            if (i < 1000 - 1) {
                printStream.print(fragSizesChimericPairs[i] + ", ");
            } else {
                printStream.print(fragSizesChimericPairs[i] + "\n");
            }
        }
        printStream.print("\n");
        printStream.print("chimeric_fragment_size_active_count_array:");
        for(int i=0; i<1000; i++) {
            if (i < 1000 - 1) {
                printStream.print(fragSizesActiveChimericPairs[i] + ", ");
            } else {
                printStream.print(fragSizesActiveChimericPairs[i] + "\n");
            }
        }
        printStream.print("\n");
        printStream.print("un_ligated_fragment_size_count_array:");
        for(int i=0; i<1000; i++) {
            if (i < 1000 - 1) {
                printStream.print(fragSizesUnLigatedPairs[i] + ", ");
            } else {
                printStream.print(fragSizesUnLigatedPairs[i] + "\n");
            }
        }
        printStream.print("\n");
        printStream.print("self_ligated_fragment_size_count_array:");
        for(int i=0; i<FRAG_SIZE_LIMIT; i++) {
            if (i < FRAG_SIZE_LIMIT - 1) {
                printStream.print(fragSizesSelfLigatedSameInternalPairs[i] + ", ");

            } else {
                printStream.print(fragSizesSelfLigatedSameInternalPairs[i]);

            }
        }
        printStream.print("\n");

        // prepare scatterplot for chromosome-wise clc against digest numbers
        printStream.print("\n");
        printStream.print("trans_cis_scatter_values_array:[");
        int cnt = 0;
        int trans_cnt=0;
        int cis_cnt=0;
        for (String chromosome : transCounts.keySet()) {
            if(chromosome.equals("chrM") || chromosome.equals("chrY")) {continue;}
            double chr_clc = 1.0*transCounts.get(chromosome)/(cisCounts.get(chromosome)+transCounts.get(chromosome));
            if(cnt==0) {
                printStream.print("{\"name\"%\"" + chromosome + "\", \"x\"%" + String.format("%.2f", chr_clc) + ",\"y\"%" + digestMap.getDigestMap().get(chromosome).getNumOfDigestsForChromosome() + "}");
                cnt++;
            } else {
                printStream.print(", {\"name\"%\"" + chromosome + "\", \"x\"%" + String.format("%.2f", chr_clc) + ",\"y\"%" + digestMap.getDigestMap().get(chromosome).getNumOfDigestsForChromosome() + "}");
            }
            trans_cnt = trans_cnt + transCounts.get(chromosome);
            cis_cnt = cis_cnt + cisCounts.get(chromosome);
            //logger.trace(chromosome + "\t" + transCounts.get(chromosome) + "\t" + cisCounts.get(chromosome) + "\t" + digestMap.getDigestMap().get(chromosome).getNumOfDigestsForChromosome());
        }
        printStream.print("]\n\n");
        double global_clc = 1.0*trans_cnt/(trans_cnt + cis_cnt);
        printStream.print("global_clc:" + String.format("%.4f", global_clc) + "\n");
    }

    /**
     * This function assembles the names of all summarize files.
     *
     * @param outputPathPrefix The prefix that will be used for output files, e.g.,  foo for foo.valid_pairs.aligned.bam
     */
    private void createOutputNames(String outputPathPrefix) {
        outputBAMvalid = String.format("%s.%s", outputPathPrefix, "valid_pairs.aligned.bam");
        outputBAMrejected = String.format("%s.%s", outputPathPrefix, "rejected_pairs.aligned.bam");
        outputFragSizesCountsRscript = String.format("%s.%s", outputPathPrefix, "frag.sizes.counts.script.R");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "align.stats.txt");
    }
}
