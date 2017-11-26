package org.jax.diachromatic.truncation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import org.jax.diachromatic.util.Pair;



/**
 * Parse paired end FASTQ read files. This class will read un g-zipped files.
 */
public class FastqPairParser {
    private static final Logger logger = LogManager.getLogger();
    private final String fastqFile1;
    private final String fastqFile2;

    private int nReadsProcessed;
    private int nReadsTruncated;


    /** A buffered reader for {@link #fastqFile1}.*/
    private BufferedReader br1;
    /** A buffered reader for {@link #fastqFile2}.*/
    private BufferedReader br2;

    private Pair<FastQRecord,FastQRecord> currentPair=null;

    private String ligationSequence=null;

    public FastqPairParser(String file1, String file2, String ligationSequence) {
        fastqFile1=file1;
        fastqFile2=file2;
        logger.trace(String.format("Processing FASTQ files %s and %s",file1,file2 ));
        logger.trace(String.format("Using ligation sequence %s",ligationSequence ));
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
        this.ligationSequence=ligSeq;
        movePairIterator();
    }

    public int getnReadsProcessed() {
        return nReadsProcessed;
    }

    public int getnReadsTruncated() {
        return nReadsTruncated;
    }

    /**
     * This function will read four lines each from {@link #br1} and {@link #br2}. If there is a problem, it will
     * set {@link #currentPair} to null and return. Otherwise, it will use the four lines from each file handle to
     * create two {@link FastQRecord} objects and place them in {@link #currentPair}. If there is any exception thrown
     * while constructing the {@link FastQRecord} objects, the function will set {@link #currentPair} to null and return
     * @throws IOException
     */
    private void movePairIterator() throws IOException {
        String read1[] = new String[4];
        String read2[] = new String[4];
        String line=null;
        if ((line=br1.readLine())!=null) {
            read1[0]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br1.readLine())!=null) {
            read1[1]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br1.readLine())!=null) {
            read1[2]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br1.readLine())!=null) {
            read1[3]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br2.readLine())!=null) {
            read2[0]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br2.readLine())!=null) {
            read2[1]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br2.readLine())!=null) {
            read2[2]=line;
        } else {
            currentPair=null;
            return;
        }
        if ((line=br2.readLine())!=null) {
            read2[3]=line;
        } else {
            currentPair=null;
            return;
        }
        // Now construct the FastQRecord objects
        FastQRecord fqr1 = new FastQRecord(read1);
        FastQRecord fqr2 = new FastQRecord(read2);
        boolean a= fqr1.truncateIfLigationSiteFound();
        boolean b =fqr2.truncateIfLigationSiteFound();
        nReadsProcessed++;
        if (a || b) nReadsTruncated++;


        currentPair = new Pair<>(fqr1,fqr2);
    }


    public boolean hasNextPair() {
        return currentPair!=null;
    }


    public Pair<FastQRecord,FastQRecord> getNextPair() {
        Pair<FastQRecord,FastQRecord> tmp = currentPair;
        try {
            movePairIterator();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return tmp;
    }




}
