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
     * Will be incremented for each query.
     */
    private Integer query_num;

    /**
     * Number of chromsome name pair keys. Should not become too large, because a good proportion of read pairs
     * maps to the same chromosome.
     */
    private Integer key_num;

    /**
     * Core structure of this class is a HashMap
     *
     * Keys: pairs of chromosome names, e.g. chr1:chr2 (lexicographical smaller comes always first)
     * Values: Another HashMap. Key is the coordinate of the lexicographical smaller reference (chr1). Value is a set
     * of integers, which are the coordinates of the lexicographical larger references (chr2).
     */
     private HashMap<String, HashMap<Integer,Set<Integer>>> dedupmap;

     DeDupMap() {
        dedupmap = new  HashMap<String, HashMap<Integer,Set<Integer>>>();
        query_num = 0;
     }

     public boolean hasSeen(ReadPair readPair) {

         query_num++;

         String stringKey;
         Integer intKey;
         Integer intVal;

         // lexicographical smaller reference ID comes always first
         if(readPair.getFivePrimeEndPosOfR1()<readPair.getFivePrimeEndPosOfR1()) {
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

         if(dedupmap.containsKey(stringKey)) {
             // a read pair has already been seen for this pair of chromosomes
             if(dedupmap.get(stringKey).containsKey(intKey)) {
                 // the coordinate of the lexicographical smaller reference ID has already been seen
                 if(dedupmap.get(stringKey).get(intKey).contains(intVal)) {
                     // the coordinate of the lexicographical smaller reference ID has already been seen
                     return true;
                 }
                 else {
                     // the coordinate of the lexicographical smaller reference ID has not been seen -> add to set
                     dedupmap.get(stringKey).get(intKey).add(intVal);
                     return false;
                 }
             }
             else {
                 // the coordinate of the lexicographical smaller reference ID has not been seen
                 Set<Integer> newSet = new HashSet<>(); // create a new integer set
                 newSet.add(intVal); // add coordinate to set
                 dedupmap.get(stringKey).put(intKey,newSet); // put set on integer HashMap
                 return false;
             }
         } else {
             // no read pair has already been seen for this pair of chromosomes
             HashMap newIntHashMap = new HashMap<Integer,Set<Integer>>(); // create new integer HashMap
             Set newSet = new HashSet<Integer>(); // create new integer set
             newSet.add(intVal); /// add coordinate to set
             newIntHashMap.put(stringKey,newSet);  // put set on integer HashMap
             dedupmap.put(stringKey,newIntHashMap); // put integer HashMap on string HashMap
             return false;
         }
     }
}