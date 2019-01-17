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
     private HashMap<String, HashMap<Integer,Set<Integer>>> dedupmap;

     DeDupMap(boolean useRelativeOrientation) {
        dedupmap = new  HashMap<String, HashMap<Integer,Set<Integer>>>();
        query_num = 0;
        insertion_num = 0;
        first_coord_num = 0;
        second_coord_num = 0;
        this.useRelativeOrientation=useRelativeOrientation;
     }

     public boolean hasSeen(ReadPair readPair) {

         query_num++;

         String stringKey;
         Integer intKey;
         Integer intVal;

         if(readPair.isTrans()) {
            // if reads are on different chromosomes, read with lexicographically smaller reference ID comes always first
            if(readPair.forward().getReferenceName().compareTo(readPair.reverse().getReferenceName())<0) {
                 stringKey = readPair.forward().getReferenceName();
                 stringKey += readPair.reverse().getReferenceName();
                 intKey = readPair.getFivePrimeEndPosOfR1();
                 intVal = readPair.getFivePrimeEndPosOfR2();
             }
             else {
                 stringKey = readPair.reverse().getReferenceName();
                 stringKey += readPair.forward().getReferenceName();
                 intKey = readPair.getFivePrimeEndPosOfR2();
                 intVal = readPair.getFivePrimeEndPosOfR1();
            }
         } else {
             // if reads are on the same chromosome, read with smaller coordinate ID comes always first
             if(readPair.getFivePrimeEndPosOfR1()<readPair.getFivePrimeEndPosOfR2()) {
                 stringKey = readPair.forward().getReferenceName();
                 stringKey += readPair.reverse().getReferenceName();
                 intKey = readPair.getFivePrimeEndPosOfR1();
                 intVal = readPair.getFivePrimeEndPosOfR2();
             }
             else {
                 stringKey = readPair.forward().getReferenceName();
                 stringKey += readPair.reverse().getReferenceName();
                 intKey = readPair.getFivePrimeEndPosOfR2();
                 intVal = readPair.getFivePrimeEndPosOfR1();
             }
         }

         /*
          * Also care about relative orientation of read pair, i.e. two read pairs with identical coordinates
          * (reference names and 5' end positions) but different relative orientation are not regarded as
          * duplicated.
          */
         if(useRelativeOrientation) {
             if (readPair.getRelativeOrientationTag().equals("F1R2") || readPair.getRelativeOrientationTag().equals("F2R1")) {
                 stringKey += "1";
             }
             if (readPair.getRelativeOrientationTag().equals("R1F2") || readPair.getRelativeOrientationTag().equals("R2F1")) {
                 stringKey += "2";
             }
             if (readPair.getRelativeOrientationTag().equals("F1F2") || readPair.getRelativeOrientationTag().equals("F2F1")) {
                 stringKey += "3";
             }
             if (readPair.getRelativeOrientationTag().equals("R1R2") || readPair.getRelativeOrientationTag().equals("R2R1")) {
                 stringKey += "4";
             }
         }

         if(dedupmap.containsKey(stringKey)) {
             // a read pair has already been seen for this pair of chromosomes
             if(dedupmap.get(stringKey).containsKey(intKey)) {
                 // the coordinate of the first read has already been seen
                 if(dedupmap.get(stringKey).get(intKey).contains(intVal)) {
                     // the coordinate of the second read has already been seen
                     return true;
                 }
                 else {
                     // coordinate of the second read has not yet been seen
                     dedupmap.get(stringKey).get(intKey).add(intVal); // add to set
                     second_coord_num++;
                     insertion_num++;
                     logger.trace(stringKey + "\t" + intKey + "\t" + intVal);
                     return false;
                 }
             }
             else {
                 // the coordinate of the first read has not yet been seen
                 Set<Integer> newSet = new HashSet<>(); // create a new integer set
                 newSet.add(intVal); // add coordinate to set
                 dedupmap.get(stringKey).put(intKey,newSet); // put set on integer HashMap
                 first_coord_num++;
                 second_coord_num++;
                 insertion_num++;
                 //logger.trace(stringKey + "\t" + intKey + "\t" + intVal);
                 return false;
             }
         } else {
             // a read pair has not yet been seen for this pair of chromosomes
             HashMap<Integer,Set<Integer>> newIntHashMap = new HashMap<>(); // create new integer HashMap
             Set newSet = new HashSet<Integer>(); // create new integer set
             newSet.add(intVal); /// add coordinate to set
             newIntHashMap.put(intKey,newSet);  // put set on integer HashMap
             dedupmap.put(stringKey,newIntHashMap); // put integer HashMap on string HashMap
             first_coord_num++;
             second_coord_num++;
             insertion_num++;
             logger.trace(stringKey + "\t" + intKey + "\t" + intVal);
             return false;
         }
     }

     /** @return Number of keys used to kep track of chromosome interaction combinations, e.g., chr1chr23 means chromosome 1-chromosome2-orientation 3.*/
     public int getNumOfChrPairKeys() {
         return dedupmap.size();
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
}