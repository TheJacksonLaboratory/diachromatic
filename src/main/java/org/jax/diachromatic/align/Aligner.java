package org.jax.diachromatic.align;


import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileWriterFactory;


import htsjdk.samtools.util.ProgressLogger;
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
 * Note that we have made several of the functions in this class package access for testing purposes
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.2 (2018-01-06)
 */
public class Aligner {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);


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
    private int n_paired = 0;
    private int n_paired_unique = 0;

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

    private int  n_unmappedPair = 0;

    private String outputBAMvalid, outputBAMrejected, outputTsvInteractingFragmentCounts, outputTsvInteractionCounts, outputFragSizesCountsRscript, outputTxtStats;

    private String filenamePrefix;

    private int n_could_not_assign_to_digest = 0;

    /**
     * Number of unique trans pairs.
     */
    private int n_trans_pairs = 0;

    /**
     * HashSet to keep track of duplicates.
     */
    private HashSet deDupHash = new HashSet();

    /** Number of read pairs whose insert was found to have a size above the threshold defined  in {@link ReadPair}.*/

    private int n_same_internal = 0;
    private int n_same_dangling_end = 0;
    private int n_same_circularized_internal = 0;
    private int n_same_circularized_dangling = 0;
    private int n_religation = 0;
    private int n_contiguous = 0;
    private int n_insert_too_short = 0;
    private int n_insert_too_long = 0;
    private int n_valid_pairs =0;
    private int n_not_categorized=0;

    /* count variables for different orientations of read pairs */

    private int n_F1F2 = 0;
    private int n_F2F1 = 0;
    private int n_R1R2 = 0;
    private int n_R2R1 = 0;
    private int n_F1R2 = 0;
    private int n_R1F2 = 0;
    private int n_R2F1 = 0;
    private int n_F2R1 = 0;

    private int n_F1F2_valid = 0;
    private int n_F2F1_valid = 0;
    private int n_R1R2_valid = 0;
    private int n_R2R1_valid = 0;
    private int n_F1R2_valid = 0;
    private int n_R1F2_valid = 0;
    private int n_R2F1_valid = 0;
    private int n_F2R1_valid = 0;

    private static int FRAG_SIZE_LIMIT = 10000;
    private int[] fragSizesAllPairs =  new int[FRAG_SIZE_LIMIT+1];
    private int[] fragSizesHybridActivePairs =  new int[FRAG_SIZE_LIMIT+1];


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
    private  DigestMap digestMap;

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

    /** Counter up the number of errors encountered in our reads. THe key is the type of error, and the value is
     * the count over the entire pair of SAM files.
     */
    private Map<ErrorCode,Integer> errorCounts;


    /**
     * Handle to write valid reads. Used in {@link #inputSAMfiles()}.
     */
    private SAMFileWriter validReadsWriter;
    /**
     * Handle to write invalid reads. Used in {@link #inputSAMfiles()}.
     */
    private SAMFileWriter rejectedReadsWriter;

    /**
     * If set to true, rejected readpairs are output to {@link #outputRejectedReads} .
     */
    private final boolean outputRejectedReads;

    private Integer lowerFragSize;
    private Integer upperFragSize;

    private boolean useStringentUniqueSettings;

    /**
     * Stores counts about interactions
     */
    InteractionCountsMap interactionMap;

    /**
     * @param sam1    SAM file for the truncated "forward" reads
     * @param sam2    SAM file for the truncated "reverse" reads
     * @param digests see {@link #digestmap}.
     */
    public Aligner(String sam1, String sam2, Map<String, List<Digest>> digests, boolean outputRejected, String outputPathPrefix, DigestMap digestMap, Integer lowerFragSize, Integer upperFragSize, String filenamePrefix, boolean useStringentUniqueSettings) {
        samPath1 = sam1;
        samPath2 = sam2;
        reader1 = SamReaderFactory.makeDefault().open(new File(samPath1));
        reader2 = SamReaderFactory.makeDefault().open(new File(samPath2));
        it1 = reader1.iterator();
        it2 = reader2.iterator();
        digestmap = digests;
        this.digestMap = digestMap;
        outputRejectedReads = outputRejected;
        this.lowerFragSize=lowerFragSize;
        this.upperFragSize=upperFragSize;
        this.filenamePrefix=filenamePrefix;
        this.useStringentUniqueSettings=useStringentUniqueSettings;

        Arrays.fill(fragSizesAllPairs, 0);
        Arrays.fill(fragSizesHybridActivePairs, 0);

        VERSION = Commandline.getVersion();
        initializeErrorMap();
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
            return new ReadPair(record1, record2, digestmap, digestMap, lowerFragSize, upperFragSize, useStringentUniqueSettings);
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

        SAMFileHeader header = reader1.getFileHeader();
        SAMFileHeader header2 = reader2.getFileHeader();
        // first add program records from the reverse SAM file
        List<SAMProgramRecord> pgList = header2.getProgramRecords();
        for (SAMProgramRecord spr : pgList) {
            //header.addProgramRecord(spr);
        }
        // now add the new program record from Diachromatic
        String programGroupId = "@PG\tID:Diachromatic\tPN:Diachromatic\tVN:" + VERSION;
        SAMProgramRecord programRecord = new SAMProgramRecord(programGroupId);
        header.addProgramRecord(programRecord);
        // we are good to go with this SAMFileHeader

        // init BAM outfile
        boolean presorted = false;
        this.validReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMvalid));
        if(outputRejectedReads) {
            this.rejectedReadsWriter = new SAMFileWriterFactory().makeBAMWriter(header, presorted, new File(outputBAMrejected));
        }
        final ProgressLogger pl = new ProgressLogger(log, 1000000);

        interactionMap = new InteractionCountsMap(1);

        DeDupMap deDupMapAll = new DeDupMap(true);

        ReadPair pair;

        while ((pair = getNextPair())!= null) {

            n_total++;

            // first check whether both reads were mapped
            if(pair.isUnMappedR1()) {n_unmapped_read1++;}
            if(pair.isUnMappedR2()) {n_unmapped_read2++;}
            if(pair.isUnMappedR1()||pair.isUnMappedR2()) {n_unmappedPair++;}
            if(pair.isMultiMappedR1()) {n_multimapped_read1++;}
            if(pair.isMultiMappedR2()) {n_multimapped_read2++;}
            if(pair.isMultiMappedR1()||pair.isMultiMappedR2()) {n_multimappedPair++;}

            // count categories of pairs
            if(pair.isPaired()) {

                n_paired++;

                // de-duplication starts with paired pairs
                if(deDupMapAll.hasSeen(pair)) {
                    n_duplicate++;
                    continue;
                }

                n_paired_unique++;

                if(pair.getCategoryTag().equals("SI")) {n_same_internal++;}
                if(pair.getCategoryTag().equals("DE")) {n_same_dangling_end++;}
                if(pair.getCategoryTag().equals("CI")) {n_same_circularized_internal++;}
                if(pair.getCategoryTag().equals("CD")) {n_same_circularized_dangling++;}
                if(pair.getCategoryTag().equals("RL")) {n_religation++;}
                if(pair.getCategoryTag().equals("CT")) {n_contiguous++;}
                if(pair.getCategoryTag().equals("TS")) {n_insert_too_short++;}
                if(pair.getCategoryTag().equals("TL")) {n_insert_too_long++;}
                if(pair.getCategoryTag().equals("VP")) {n_valid_pairs++;}
                if(pair.getCategoryTag().equals("NA")) {n_not_categorized++;}

                if(pair.isTrans()) {
                    n_trans_pairs++;
                }
            }

            // both reads were uniquely mapped, otherwise continue
            if(!pair.isPaired()) {updateErrorMap(pair.getErrorCodes()); continue;}


            if(pair.getRelativeOrientationTag().equals("F1F2")) {n_F1F2++;}
            if(pair.getRelativeOrientationTag().equals("F2F1")) {n_F2F1++;}
            if(pair.getRelativeOrientationTag().equals("R1R2")) {n_R1R2++;}
            if(pair.getRelativeOrientationTag().equals("R2R1")) {n_R2R1++;}
            if(pair.getRelativeOrientationTag().equals("F1R2")) {n_F1R2++;}
            if(pair.getRelativeOrientationTag().equals("R2F1")) {n_R2F1++;}
            if(pair.getRelativeOrientationTag().equals("F2R1")) {n_F2R1++;}
            if(pair.getRelativeOrientationTag().equals("R1F2")) {n_R1F2++;}

            if(pair.isValid()) {
                if(pair.getRelativeOrientationTag().equals("F1F2")) {n_F1F2_valid++;}
                if(pair.getRelativeOrientationTag().equals("F2F1")) {n_F2F1_valid++;}
                if(pair.getRelativeOrientationTag().equals("R1R2")) {n_R1R2_valid++;}
                if(pair.getRelativeOrientationTag().equals("R2R1")) {n_R2R1_valid++;}
                if(pair.getRelativeOrientationTag().equals("F1R2")) {n_F1R2_valid++;}
                if(pair.getRelativeOrientationTag().equals("R2F1")) {n_R2F1_valid++;}
                if(pair.getRelativeOrientationTag().equals("F2R1")) {n_F2R1_valid++;}
                if(pair.getRelativeOrientationTag().equals("R1F2")) {n_R1F2_valid++;}
            }


            // count sizes of all fragments
            Integer incrementFragSize = pair.getCalculatedInsertSize();
            if(FRAG_SIZE_LIMIT<incrementFragSize) { incrementFragSize = FRAG_SIZE_LIMIT; }
            fragSizesAllPairs[incrementFragSize]++;

            // count sizes of all hybrid active fragments
            if((pair.forwardDigestIsActive() & !pair.reverseDigestIsActive()) || (!pair.forwardDigestIsActive() & pair.reverseDigestIsActive())) {
                fragSizesHybridActivePairs[incrementFragSize]++;
            }

            if(pair.isValid()){

                validReadsWriter.addAlignment(pair.forward());
                validReadsWriter.addAlignment(pair.reverse());

                n_good++;

                // count interaction
                /*interactionMap.incrementFragPair(0,
                        pair.forward().getReferenceName(),
                        pair.getForwardDigestStart(),
                        pair.getForwardDigestEnd(),
                        pair.forwardDigestIsActive(),
                        pair.reverse().getReferenceName(),
                        pair.getReverseDigestStart(),
                        pair.getReverseDigestEnd(),
                        pair.reverseDigestIsActive(),
                        pair.getRelativeOrientationTag());*/
                Integer relOriTag=-1;
                if(pair.getRelativeOrientationTag().equals("F1F2")) {relOriTag=0;} // twisted
                if(pair.getRelativeOrientationTag().equals("F2F1")) {relOriTag=1;} // twisted
                if(pair.getRelativeOrientationTag().equals("R1R2")) {relOriTag=2;} // twisted
                if(pair.getRelativeOrientationTag().equals("R2R1")) {relOriTag=3;} // twisted
                if(pair.getRelativeOrientationTag().equals("F1R2")) {relOriTag=4;} // simple
                if(pair.getRelativeOrientationTag().equals("R2F1")) {relOriTag=5;} // simple
                if(pair.getRelativeOrientationTag().equals("F2R1")) {relOriTag=6;} // simple
                if(pair.getRelativeOrientationTag().equals("R1F2")) {relOriTag=7;} // simple
                interactionMap.incrementFragPair2(0, pair.getForwardDigestKey(), pair.getReverseDigestKey(),relOriTag);
                if(interactionMap.getTotalNumberOfInteractionsForCondition(0)%10000==0) { logger.trace("Size of interactionMap: " + interactionMap.getTotalNumberOfInteractionsForCondition(0)); }
            } else {
                updateErrorMap(pair.getErrorCodes());
                if (outputRejectedReads) {
                    rejectedReadsWriter.addAlignment(pair.forward());
                    rejectedReadsWriter.addAlignment(pair.reverse());
                }
                // discard this read and go to the next one
            }
        }
        logger.trace(outputTsvInteractionCounts);
        //interactionMap.printInteractionCountsMapAsCountTable(outputTsvInteractionCounts);
        interactionMap.printInteractionCountsMapAsCountTable2(digestMap,outputTsvInteractionCounts);

        //interactionMap.printFragmentInteractionCountsMapAsCountTable(outputTsvInteractingFragmentCounts);
        validReadsWriter.close();
        if(outputRejectedReads) {
            rejectedReadsWriter.close();
        }

        printFragmentLengthDistributionRscript(fragSizesAllPairs, fragSizesHybridActivePairs);

        logger.trace("" );
        logger.trace("Deduplication stats:" );
        logger.trace("n_duplicate: " + n_duplicate);
        logger.trace("deDupMapAll.getNumOfChrPairKeys(): " + deDupMapAll.getNumOfChrPairKeys());
        logger.trace("deDupMapAll.getNumOfQueries(): " + deDupMapAll.getNumOfQueries());
        logger.trace("deDupMapAll.getNumOfInsertions(): " + deDupMapAll.getNumOfInsertions());
        logger.trace("deDupMapAll.getNumOfFirstCoords(): " + deDupMapAll.getNumOfFirstCoords());
        logger.trace("deDupMapAll.getNumOfSecondCoords(): " + deDupMapAll.getNumOfSecondCoords());
        logger.trace("" );







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

        printStream.print("plot(length,fragSizesHybridActivePairs,main=MAIN, xlim=XLIM,type=\"l\", ylim=c(0,YLIM),col=\"red\",ylab=\"fragment count\")\n");

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

      /*
        MEAN_FRAG_SIZE_ALL<-sum(length*fragSizesAllPairs)/sum(fragSizesAllPairs)
        print(MEAN_FRAG_SIZE_ALL)

        MEAN_FRAG_SIZE_ACTIVE<-sum(length*fragSizesHybridActivePairs)/sum(fragSizesHybridActivePairs)
        print(MEAN_FRAG_SIZE_ACTIVE)

        s<-0
        MEADIAN_FRAG_SIZE<-1
        while(s < sum(fragSizesAllPairs)/2) {
            s = s + fragSizesAllPairs[MEADIAN_FRAG_SIZE]
            MEADIAN_FRAG_SIZE = MEADIAN_FRAG_SIZE + 1
        }
        print(MEADIAN_FRAG_SIZE)

        s<-0
        MEADIAN_ACTIVE_FRAG_SIZE<-1
        while(s < sum(fragSizesAllPairs)/2) {
            s = s + fragSizesAllPairs[MEADIAN_ACTIVE_FRAG_SIZE]
            MEADIAN_ACTIVE_FRAG_SIZE = MEADIAN_ACTIVE_FRAG_SIZE + 1
        }
        print(MEADIAN_ACTIVE_FRAG_SIZE)
    */
    }



    /** The align {@link #errorCounts} is initialize by setting the counts for all elements to zero. */
    private void initializeErrorMap() {
        this.errorCounts=new HashMap<>();
        for (ErrorCode ec : ErrorCode.values()) {
            errorCounts.put(ec,0);
        }
    }

    /**
     * Increment the error code for the errors encounted in a read pair.
     * @param errors Set of errors encountered for some read pair.
     */
    private void updateErrorMap(Set<ErrorCode> errors) {
        for (ErrorCode ec : errors) {
            errorCounts.put(ec,1+errorCounts.get(ec));
        }
    }


    public void printStatistics() throws FileNotFoundException {

        logger.trace(String.format("n_total pairs=%d\n", n_total));

        logger.trace(String.format("n_unmapped_read1=%d", n_unmapped_read1));
        logger.trace(String.format("n_unmapped_read2=%d\n", n_unmapped_read2));

        logger.trace(String.format("n_multimapped_read1=%d", n_multimapped_read1));
        logger.trace(String.format("n_multimapped_read2=%d\n", n_multimapped_read2));
        logger.trace(String.format("n_multimappedPair=%d\n", n_multimappedPair));


        logger.trace(String.format("n_paired=%d (%.1f%%)\n", n_paired, (100.0 * n_paired / n_total)));

        logger.trace(String.format("n_could_not_assign_to_digest=%d (%.1f%%)\n", n_could_not_assign_to_digest, (100.0 * n_could_not_assign_to_digest / n_total)));

        logger.trace(String.format("n_same_internal=%d", n_same_internal));
        logger.trace(String.format("n_same_dangling_end=%d", n_same_dangling_end));
        logger.trace(String.format("n_same_circularized_read=%d", n_same_circularized_internal+n_same_circularized_dangling));
        logger.trace(String.format("n_religation=%d", n_religation));
        logger.trace(String.format("n_contiguous=%d\n", n_contiguous));
        logger.trace(String.format("n_not_categorized=%d\n", n_not_categorized));

        logger.trace(String.format("n_insert_too_long=%d  (%.1f%%)", n_insert_too_long, (100.0 * n_insert_too_long / n_paired_unique)));
        logger.trace(String.format("n_insert_too_short=%d  (%.1f%%)\n", n_insert_too_short, (100.0 * n_insert_too_short / n_paired_unique)));

        logger.trace(String.format("n_valid_pairs=%d (%.1f%%)", n_valid_pairs, (100.0 * n_valid_pairs / n_paired_unique)));
        logger.trace("");
        logger.trace("Total number of pairs: " + (n_same_internal+n_same_dangling_end+n_same_circularized_internal+n_same_circularized_dangling+n_religation+n_contiguous+n_insert_too_long+n_insert_too_short+n_valid_pairs));

        logger.trace("");
        logger.trace("Distribution of pair orientations (all pairs):");
        logger.trace("n_F1F2 (commie)" + "\t" + n_F1F2);
        logger.trace("n_F2F1 (commie)" + "\t" +  n_F2F1);
        logger.trace("n_R1R2 (commie)" + "\t" +  n_R1R2);
        logger.trace("n_R2R1 (commie)" + "\t" +  n_R2R1);
        logger.trace("n_F1R2 (innie)" + "\t" +  n_F1R2);
        logger.trace("n_F2R1 (innie)" + "\t" +  n_F2R1);
        logger.trace("n_R2F1 (outie)" + "\t" +  n_R2F1);
        logger.trace("n_R1F2 (outie)" + "\t" +  n_R1F2);
        logger.trace("");
        Integer n_outies = n_R1F2 + n_R2F1;
        Integer n_innies = n_F1R2 + n_F2R1;
        Integer n_commies = n_F1F2 + n_R1R2 + n_R2R1 + n_F2F1;
        logger.trace("n_outies: " + n_outies);
        logger.trace("n_innies: " + n_innies);
        logger.trace("n_commies: " + n_commies);


        /**
         * Self-ligation must only result in read pairs that are pointing in outward direction (outies).
         * Simple loops can also result in (outies), whereas twisted loops exclusively result in read pairs
         * pointing in the same direction. Twisted loops are less affected by biases regarding proximity, e.g.
         * same internal artifacts are almost always innies. Therefore, the number of outies that resulted from
         * self-ligation is estimated by subtracting the mean number of commies from the mean number of commies.
         * The proportion of outies that resulted from self-ligation and not simple loops is the self-ligation
         * coefficient. Pseudo counts were added to avoid divisions by zero in extreme cases.
         * The coefficient approaches zero, if no self-ligation occurred.
         *
         */
        float mean_outies = n_outies/2;
        float mean_commies = n_commies/4;
        logger.trace("mean_outies: " + mean_outies);
        logger.trace("mean_commies: " + mean_commies);
        float selfLigationCoefficient = (float) ((1.0*(mean_outies-mean_commies)/(mean_outies+1)));

        logger.trace("");
        logger.trace("Summary statistics about interactions between active and inactive fragments:");
        logger.trace("");
        logger.trace("\t" + "Total number of interactions: " + interactionMap.getTotalNumberOfInteractionsForCondition(0));
        logger.trace("\t" + "Number of interactions between active fragments: " + interactionMap.getNumberOfInteractionsBetweenActiveFragmentsForCondition(0));
        logger.trace("\t" + "Number of interactions between inactive fragments: " + interactionMap.getNumberOfInteractionsBetweenInactiveFragmentsForCondition(0));
        logger.trace("\t" + "Number of interactions between active and inactive fragments: " + interactionMap.getNumberOfInteractionsBetweenActiveAndInactiveFragmentsForCondition(0));
        logger.trace("");
        logger.trace("\t" + "Total number of interacting fragments: " + interactionMap.getTotalNumberOfInteractingFragmentsForCondition(0));
        logger.trace("\t" + "Number of active interacting fragments: " + interactionMap.getTotalNumberOfActiveInteractingFragmentsForCondition(0));
        logger.trace("");
        logger.trace("\t" + "Enrichment Coefficients:");
        logger.trace("\t\t" + "Target Enrichment Coefficient (TEC): " + String.format("%.2f%%", 100*interactionMap.getTargetEnrichmentCoefficientForCondition(0)));
        logger.trace("\t\t" + "Valid Interaction Enrichment Coefficient (VIEC): " + String.format("%.2f%%", 100.0*n_valid_pairs/n_total));
        logger.trace("\t\t" + "Self Ligation Coefficient (SLC): " + String.format("%.2f", selfLigationCoefficient));
        logger.trace("\t\t" + "Cross-ligation coefficient (CLC): " + String.format("%.2f%%", 100.0*n_trans_pairs/n_paired_unique));
        logger.trace("\t\t" + "Re-ligation coefficient (RLC): " + String.format("%.2f%%", 100.0*(n_paired_unique-n_same_dangling_end)/n_paired_unique));
        logger.trace("\t\t" + "Pair duplication rate: " + String.format("%.2f%%", 100.0*n_duplicate/n_paired));
        logger.trace("n_duplicate: " + n_duplicate);
        logger.trace("");

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(outputTxtStats));

        printStream.print("Summary statistics\n");
        printStream.print("==================\n\n");
        printStream.print("\n");
        printStream.print("Alignment statistics\n");
        printStream.print("--------------------\n");
        printStream.print("\n");
        printStream.print("Total number of read pairs processed:\t" + n_total + "\n");

        printStream.print("Number of unmapped read pairs:\t" + n_unmappedPair + String.format(" (%.2f%%)", 100.0*n_unmappedPair/n_total) + "\n");
        printStream.print("\tNumber of unpapped R1 reads:\t" + n_unmapped_read1 + "\n");
        printStream.print("\tNumber of unpapped R2 reads:\t" + n_unmapped_read2 + "\n");

        printStream.print("Number of multimapped read pairs:\t" + n_multimappedPair + String.format(" (%.2f%%)", 100.0*n_multimappedPair/n_total) + "\n");
        printStream.print("\tNumber of multimapped R1 reads:\t" + n_multimapped_read1 + "\n");
        printStream.print("\tNumber of multimapped R2 reads:\t" + n_multimapped_read2 + "\n");
        printStream.print("Note:\tThere may an overlap between unmapped and multimapped pairs." + "\n");

        printStream.print("Number of paired read pairs:\t" + n_paired + String.format(" (%.2f%%)", 100.0*n_paired/n_total) + "\n");
        printStream.print("\tNumber of unique paired read pairs:\t" + n_paired_unique + "\n");
        printStream.print("\tNumber of duplicated pairs:\t" + n_duplicate + "\n");

        printStream.print("\tPair duplication rate:\t" + String.format("%.2f%%", 100.0*n_duplicate/n_paired) + "\n");
        printStream.print("\n");
        printStream.print("Read pair orientations of unique paired read pairs:\n");
        printStream.print("F1F2 - commie:\t" + n_F1F2 + String.format(" (%.2f%%)", 100.0*n_F1F2/n_paired_unique) + "\n");
        printStream.print("F2F1 - commie:\t" + n_F2F1 + String.format(" (%.2f%%)", 100.0*n_F2F1/n_paired_unique) + "\n");
        printStream.print("R1R2 - commie:\t" + n_R1R2 + String.format(" (%.2f%%)", 100.0*n_R1R2/n_paired_unique) + "\n");
        printStream.print("R2R1 - commie:\t" + n_R2R1 + String.format(" (%.2f%%)", 100.0*n_R2R1/n_paired_unique) + "\n");
        printStream.print("F1R2 - innie:\t" + n_F1R2 + String.format(" (%.2f%%)", 100.0*n_F1R2/n_paired_unique) + "\n");
        printStream.print("F2R1 - innie:\t" + n_F2R1 + String.format(" (%.2f%%)", 100.0*n_F2R1/n_paired_unique) + "\n");
        printStream.print("R2F1 - outie:\t" + n_R2F1 + String.format(" (%.2f%%)", 100.0*n_R2F1/n_paired_unique) + "\n");
        printStream.print("R1F2 - outie:\t" + n_R1F2 + String.format(" (%.2f%%)", 100.0*n_R1F2/n_paired_unique) + "\n");
        printStream.print("\n");
        printStream.print("Read pair orientations of unique valid paired read pairs:\n");
        printStream.print("F1F2 - commie:\t" + n_F1F2_valid + String.format(" (%.2f%%)", 100.0*n_F1F2_valid/n_valid_pairs) + "\n");
        printStream.print("F2F1 - commie:\t" + n_F2F1_valid + String.format(" (%.2f%%)", 100.0*n_F2F1_valid/n_valid_pairs) + "\n");
        printStream.print("R1R2 - commie:\t" + n_R1R2_valid + String.format(" (%.2f%%)", 100.0*n_R1R2_valid/n_valid_pairs) + "\n");
        printStream.print("R2R1 - commie:\t" + n_R2R1_valid + String.format(" (%.2f%%)", 100.0*n_R2R1_valid/n_valid_pairs) + "\n");
        printStream.print("F1R2 - innie:\t" + n_F1R2_valid + String.format(" (%.2f%%)", 100.0*n_F1R2_valid/n_valid_pairs) + "\n");
        printStream.print("F2R1 - innie:\t" + n_F2R1_valid + String.format(" (%.2f%%)", 100.0*n_F2R1_valid/n_valid_pairs) + "\n");
        printStream.print("R2F1 - outie:\t" + n_R2F1_valid + String.format(" (%.2f%%)", 100.0*n_R2F1_valid/n_valid_pairs) + "\n");
        printStream.print("R1F2 - outie:\t" + n_R1F2_valid + String.format(" (%.2f%%)", 100.0*n_R1F2_valid/n_valid_pairs) + "\n");
        printStream.print("\n");
        printStream.print("\n");
        printStream.print("Artifact statistics (HiCUP like)\n");
        printStream.print("--------------------------------\n");
        printStream.print("\n");
        printStream.print("Valid:\t" + n_valid_pairs + String.format(" (%.2f%%)", 100.0*n_valid_pairs/n_paired_unique) + "\n");
        printStream.print("Same internal:\t" + n_same_internal + String.format(" (%.2f%%)", 100.0*n_same_internal/n_paired_unique) + "\n");
        printStream.print("Same dangling:\t" + n_same_dangling_end + String.format(" (%.2f%%)", 100.0*n_same_dangling_end/n_paired_unique) + "\n");
        Integer n_same_circularized = n_same_circularized_internal+n_same_circularized_dangling;
        printStream.print("Same circularized:\t" + n_same_circularized + String.format(" (%.2f%%)", 100.0*n_same_circularized/n_paired_unique) + "\n");
        printStream.print("Re-ligation:\t" + n_religation + String.format(" (%.2f%%)", 100.0*n_religation/n_paired_unique) + "\n");
        printStream.print("Contiguous:\t" + n_contiguous + String.format(" (%.2f%%)", 100.0*n_contiguous/n_paired_unique) + "\n");
        printStream.print("Insert too short:\t" + n_insert_too_short + String.format(" (%.2f%%)", 100.0*n_insert_too_short/n_paired_unique) + "\n");
        printStream.print("Insert too long:\t" + n_insert_too_long + String.format(" (%.2f%%)", 100.0*n_insert_too_long/n_paired_unique) + "\n");
        printStream.print("Not categorized:\t" + n_not_categorized + String.format(" (%.2f%%)", 100.0*n_not_categorized/n_paired_unique) + "\n");
        printStream.print("\n");
        printStream.print("Artifact statistics (Diachromatic like)\n");
        printStream.print("---------------------------------------\n");
        printStream.print("\n");
        printStream.print("\n");
        printStream.print("Will be added soon.\n");
        printStream.print("\n");
        printStream.print("\n");
        printStream.print("Quality metrics for experimental trouble shooting \n");
        printStream.print("--------------------------------------------------\n");
        printStream.print("\n");
        printStream.print("Yield of Valid Pairs (YVP):\t" + String.format("%.2f%%", 100.0*n_valid_pairs/n_total) + "\n");
        printStream.print("Target Enrichment Coefficient (TEC):\t" + String.format("%.2f%%", 100*interactionMap.getTargetEnrichmentCoefficientForCondition(0)) + "\n");
        printStream.print("Cross-ligation coefficient (CLC):\t" + String.format("%.2f%%", 100.0*n_trans_pairs/n_paired_unique) + "\n");
        printStream.print("Re-ligation coefficient (RLC):\t" + String.format("%.2f%%", 100.0*(n_paired_unique-n_same_dangling_end)/n_paired_unique) + "\n");
        printStream.print("\n");














    }

    private void createOutputNames(String outputPathPrefix) {
        outputBAMvalid = String.format("%s.%s", outputPathPrefix, "valid_pairs.aligned.bam");
        outputBAMrejected = String.format("%s.%s", outputPathPrefix, "rejected_pairs.aligned.bam");
        outputTsvInteractingFragmentCounts = String.format("%s.%s", outputPathPrefix, "interacting.fragments.counts.table.tsv"); // will be moved to class counts
        outputTsvInteractionCounts = String.format("%s.%s", outputPathPrefix, "interaction.counts.table.tsv"); // will be moved to class counts
        outputFragSizesCountsRscript = String.format("%s.%s", outputPathPrefix, "frag.sizes.counts.script.R");
        outputTxtStats = String.format("%s.%s", outputPathPrefix, "align.stats.txt");
    }
}
