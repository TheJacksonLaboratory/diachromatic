package org.jax.diachromatic.align;

import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class DeDupMap {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /**
     * The total number of queries to this object (performed via the {@link #hasSeen(ReadPair)} method)..
     */
    private int query_num;

    /**
     * Number of insertions that were made. Corresponds to the number of unique pairs.
     */
    private int insertion_num;

    /**
     * Number of first and second coordinates.
     */
    private int first_coord_num, second_coord_num;

    /**
     *
     */
    private boolean useRelativeOrientation;

    /**
     * Core structure of this class is a HashMap
     *
     * Keys: pairs of chromosome names, e.g. chr1:chr2 (lexicographical smaller comes always first)
     * Values: Another HashMap. Key is the coordinate of the lexicographical smaller reference (chr1). Value is a set
     * of integers, which are the coordinates of the lexicographical larger references (chr2).
     */
    private HashMap<String, Set<ReadPairCoordinates>> dedupmap2;

    DeDupMap(boolean useRelativeOrientation) {
        dedupmap2 = new  HashMap<String, Set<ReadPairCoordinates>>();
        query_num = 0;
        insertion_num = 0;
        first_coord_num = 0;
        second_coord_num = 0;
        this.useRelativeOrientation=useRelativeOrientation;
    }

    public boolean hasSeen(ReadPair readPair) {

        query_num++;

        String stringKey;
        Integer fp1;
        Integer fp2;

        if(readPair.isTrans()) {
            // if reads are on different chromosomes, read with lexicographically smaller reference ID comes always first
            if(readPair.forward().getReferenceName().compareTo(readPair.reverse().getReferenceName())<0) {
                stringKey = readPair.getReferenceSequenceOfR1();
                stringKey += readPair.getReferenceSequenceOfR2();
                fp1 = readPair.getFivePrimeEndPosOfR1();
                fp2 = readPair.getFivePrimeEndPosOfR2();

            }
            else {
                stringKey = readPair.getReferenceSequenceOfR2();
                stringKey += readPair.getReferenceSequenceOfR1();
                fp1 = readPair.getFivePrimeEndPosOfR2();
                fp2 = readPair.getFivePrimeEndPosOfR1();
            }
        } else {
            // if reads are on the same chromosome, read with smaller coordinate ID comes always first
            if(readPair.getFivePrimeEndPosOfR1()<readPair.getFivePrimeEndPosOfR2()) {
                stringKey = readPair.getReferenceSequenceOfR1();
                stringKey += readPair.getReferenceSequenceOfR2();
                fp1 = readPair.getFivePrimeEndPosOfR1();
                fp2 = readPair.getFivePrimeEndPosOfR2();
            }
            else {
                stringKey = readPair.getReferenceSequenceOfR1();
                stringKey += readPair.getReferenceSequenceOfR2();
                fp1 = readPair.getFivePrimeEndPosOfR2();
                fp2 = readPair.getFivePrimeEndPosOfR1();
            }
        }

         /*
          * Also care about relative orientation of read pair, i.e. two read pairs with identical coordinates
          * (reference names and 5' end positions) but different relative orientation are not regarded as
          * duplicated.
          */
         Integer readPairOrientation = 0;
        if(useRelativeOrientation) {
            if (readPair.getRelativeOrientationTag().equals("F1R2") || readPair.getRelativeOrientationTag().equals("F2R1")) {
                readPairOrientation = 1;
            }
            if (readPair.getRelativeOrientationTag().equals("R1F2") || readPair.getRelativeOrientationTag().equals("R2F1")) {
                readPairOrientation = 2;
            }
            if (readPair.getRelativeOrientationTag().equals("F1F2") || readPair.getRelativeOrientationTag().equals("F2F1")) {
                readPairOrientation = 3;
            }
            if (readPair.getRelativeOrientationTag().equals("R1R2") || readPair.getRelativeOrientationTag().equals("R2R1")) {
                readPairOrientation = 4;
            }
        }
        ReadPairCoordinates rpc = new ReadPairCoordinates(fp1, fp2, readPairOrientation);

        if(dedupmap2.containsKey(stringKey)) {
            // a read pair with this pair of chromosomes has already been seen
            if (dedupmap2.get(stringKey).add(rpc)) {
                insertion_num++;
                return false;
            } else {
                return true;
            }
        } else {
            // this is the first read pair with this pair of chromosomes
            Set<ReadPairCoordinates> newSet = new HashSet<>(); // create new integer set
            newSet.add(rpc);
            dedupmap2.put(stringKey,newSet);
            insertion_num++;
            return false;
        }
    }

    /** @return Number of keys used to kep track of chromosome interaction combinations, e.g., chr1chr23 means chromosome 1-chromosome2-orientation 3.*/
    public int getNumOfChrPairKeys() {
        return dedupmap2.size();
    }

    public int getNumOfQueries() {
        return query_num;
    }

    public int getNumOfInsertions() {
        return insertion_num;
    }

    public int getNumOfFirstCoords() {
        return first_coord_num;
    }

    public int getNumOfSecondCoords() {
        return second_coord_num;
    }

    public void printDeDupStatistics(int n_paired_duplicated) {
        logger.trace("" );
        logger.trace("Deduplication statitics:" );
        logger.trace("n_duplicate: " + n_paired_duplicated);
        logger.trace("getNumOfChrPairKeys(): " + getNumOfChrPairKeys());
        logger.trace("getNumOfQueries(): " + getNumOfQueries());
        logger.trace("getNumOfInsertions(): " + getNumOfInsertions());
        logger.trace("getNumOfFirstCoords(): " + getNumOfFirstCoords());
        logger.trace("getNumOfSecondCoords(): " + getNumOfSecondCoords());
        logger.trace("" );
    }
}
