package org.jax.diachromatic.truncation;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.IOException;

import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.util.Pair;


/**
 * Parse paired end FASTQ read files. This class will read g-zipped files.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
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
    /** Number of reads from file 2 that got truncated during processing. */
    private int nReadTwoTruncated;
    /** Total number of reads from each file that got processed.*/
    private int nReadsProcessed;
    /** FASTQ reader for the forward reads. */
    private FastqReader fastQreader1;
    /** FASTQ reader for the reverse reads. */
    private FastqReader fastQreader2;
    /** Number at which we show logger trace messages (every BLOCKSIZE reads)*/
    private final int BLOCKSIZE=200_000;


    /**
     * This is used to store the last processed pair of reads roughly in the style of an iterator.
     */
    private Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> currentPair = null;
    /**
     * This is the Hi-C ligation sequence that is created from the restriction enzyme.
     */
    private String ligationSequence = null;

    public FastqPairParser(String file1, String file2, String ligationSequence) throws DiachromaticException {
        fastqFile1 = file1;
        fastqFile2 = file2;
        logger.trace(String.format("Processing FASTQ files %s and %s with ligation sequence %s", file1, file2, ligationSequence));
        try {
            setUpIterator(ligationSequence);
        } catch (IOException e) {
            String msg=String.format("I/O error setting up Fastq readers for %s and %s: %s",file1,file2,e.getMessage() );
            e.printStackTrace();
            throw new DiachromaticException(msg);
        }
    }


    private void setUpIterator(String ligSeq) throws IOException {
        fastQreader1 = new FastqReader(new File(fastqFile1));
        fastQreader2 = new FastqReader(new File(fastqFile2));
        this.ligationSequence = ligSeq;
        movePairIterator();
    }

    int getnReadsProcessed() {
        return nReadsProcessed;
    }

    int getReadOneTruncated() {
        return nReadOneTruncated;
    }

    int getReadTwoTruncated() {
        return nReadTwoTruncated;
    }

    /**
     * This function will read the next FASTQ record each from {@link #fastQreader1} and {@link #fastQreader2}. If there is a problem, it will
     * set {@link #currentPair} to null and return. Otherwise, it will use the four lines from each file handle to
     * create two {@link PotentiallyTruncatedFastQRecord} objects and place them in {@link #currentPair}. If there is any exception thrown
     * while constructing the {@link PotentiallyTruncatedFastQRecord} objects, the function will set {@link #currentPair} to null and return
     *
     */
    private void movePairIterator() {
        if (fastQreader1.hasNext() && fastQreader2.hasNext()) {
            FastqRecord fq1 = fastQreader1.next();
            FastqRecord fq2 = fastQreader2.next();


            PotentiallyTruncatedFastQRecord fqr1 = new PotentiallyTruncatedFastQRecord(fq1);
            PotentiallyTruncatedFastQRecord fqr2 = new PotentiallyTruncatedFastQRecord(fq2);
            if (fqr1.truncateIfLigationSiteFound()) {
                nReadOneTruncated++;
            }
            if (fqr2.truncateIfLigationSiteFound()) {
                nReadTwoTruncated++;
            }
            currentPair = new Pair<>(fqr1, fqr2);
            nReadsProcessed++;
            if (nReadsProcessed % BLOCKSIZE == 0) {
                logger.trace("Processed readpair number {}.", nReadsProcessed);
            }
        } else {
            currentPair=null;
        }
    }

    /**
     *
     * @return true iff the FASTQ Readers have another read pair.
     */
    boolean hasNextPair() {
        return currentPair != null;
    }

    /**
     * @return the next forward/reverse read pair
     */
    Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> getNextPair() {
        Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> tmp = currentPair;
        movePairIterator();
        return tmp;
    }


}
