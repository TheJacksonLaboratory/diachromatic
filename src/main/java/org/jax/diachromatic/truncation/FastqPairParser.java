package org.jax.diachromatic.truncation;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.util.Pair;


/**
 * Parse paired end FASTQ read files. This class will read un g-zipped files.
 */
public class FastqPairParser {
    private static final Logger logger = LogManager.getLogger();
    /**
     * First read of paired end sequencing experiment.
     */
    private final String fastqFile1;
    /**
     * Second read of paired end sequencing experiment.
     */
    private final String fastqFile2;
    /**
     * Number of reads from file 1 that got truncated during processing.
     */
    private int nReadOneTruncated;
    /**
     * Number of reads from file 2 that got truncated during processing.
     */
    private int nReadTwoTruncated;
    /**
     * Total number of reads from each file that got processed.
     */
    private int nReadsProcessed;
    /**
     * A buffered reader for {@link #fastqFile1}.
     */
    private BufferedReader br1;
    /**
     * A buffered reader for {@link #fastqFile2}.
     */
    private BufferedReader br2;


    private FastqReader fastQreader1;

    private FastqReader fastQreader2;


    /**
     * This is used to store the last processed pair of reads roughly in the style of an iterator.
     */
    private Pair<FastQRecord, FastQRecord> currentPair = null;
    /**
     * This is the Hi-C ligation sequence that is created from the restriction enzyme.
     */
    private String ligationSequence = null;

    public FastqPairParser(String file1, String file2, String ligationSequence) {
        fastqFile1 = file1;
        fastqFile2 = file2;
        logger.trace(String.format("Processing FASTQ files %s and %s with ligation sequence %s", file1, file2, ligationSequence));
        try {
            setUpIterator(ligationSequence);
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1); // todo exception
        }
    }


    private void setUpIterator(String ligSeq) throws IOException {
        br1 = new BufferedReader(new FileReader(this.fastqFile1));
        br2 = new BufferedReader(new FileReader(this.fastqFile2));
//        fastQreader1 = new FastqReader(new File(fastqFile1));
//        fastQreader2 = new FastqReader(new File(fastqFile2));
        this.ligationSequence = ligSeq;
        movePairIterator();
    }

    public int getnReadsProcessed() {
        return nReadsProcessed;
    }

    public int getReadOneTruncated() {
        return nReadOneTruncated;
    }

    public int getReadTwoTruncated() {
        return nReadTwoTruncated;
    }

    /**
     * This function will read four lines each from {@link #br1} and {@link #br2}. If there is a problem, it will
     * set {@link #currentPair} to null and return. Otherwise, it will use the four lines from each file handle to
     * create two {@link FastQRecord} objects and place them in {@link #currentPair}. If there is any exception thrown
     * while constructing the {@link FastQRecord} objects, the function will set {@link #currentPair} to null and return
     *
     * @throws IOException
     */
    private void movePairIterator() throws IOException {

//        if (fastQreader2.hasNext()) {
//            FastqRecord fastqRecord1 = fastQreader2.next();
//        }


        String read1[] = new String[4];
        String read2[] = new String[4];
        String line = null;
        if (((line = br1.readLine()) != null) && line.startsWith("@")) { // TODO Error / exception if '@' not here
            read1[0] = line;
        } else {
            currentPair = null;
            return;
        }
        if ((line = br1.readLine()) != null) {
            read1[1] = line;
        } else {
            currentPair = null;
            return;
        }
        if (((line = br1.readLine()) != null) && line.startsWith("+")) {
            read1[2] = line;
        } else {
            currentPair = null;
            return;
        }
        if ((line = br1.readLine()) != null) {
            read1[3] = line;
        } else {
            currentPair = null;
            return;
        }
        if (((line = br2.readLine()) != null) && line.startsWith("@")) {
            read2[0] = line;
        } else {
            currentPair = null;
            return;
        }
        if ((line = br2.readLine()) != null) {
            read2[1] = line;
        } else {
            currentPair = null;
            return;
        }
        if (((line = br2.readLine()) != null) && line.startsWith("+")) {
            read2[2] = line;
        } else {
            currentPair = null;
            return;
        }
        if ((line = br2.readLine()) != null) {
            read2[3] = line;
        } else {
            currentPair = null;
            return;
        }
        // Now construct the FastQRecord objects
        try {
            FastQRecord fqr1 = new FastQRecord(read1);
            FastQRecord fqr2 = new FastQRecord(read2);
            if (fqr1.truncateIfLigationSiteFound()) {
                nReadOneTruncated++;
            }
            if (fqr2.truncateIfLigationSiteFound()) {
                nReadTwoTruncated++;
            }
            currentPair = new Pair<>(fqr1, fqr2);
            nReadsProcessed++;
        } catch (DiachromaticException e) {
            logger.error("Exception encountered while parsing FASTQ files: " + e.getMessage());
            e.printStackTrace();
            logger.error("Exiting program -- check the format of the input file");
            System.exit(1);
        }



    }


    public boolean hasNextPair() {
        return currentPair != null;
    }


    public Pair<FastQRecord, FastQRecord> getNextPair() {
        Pair<FastQRecord, FastQRecord> tmp = currentPair;
        try {
            movePairIterator();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return tmp;
    }


}
