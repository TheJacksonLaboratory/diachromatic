package org.jax.diachromatic.align;


import org.jax.diachromatic.exception.DiachromaticException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * This class provides a data structure that can be used get the digest coordinates associated
 * to a given coordinate x in O(log n) time.
 * <p>
 *     As input it expects to receive the digest file that was created by GOPHER. This file has the
 *     following columns
 *     <ol>
 *         <li>Chromosome, e.g., chr1</li>
 *         <li>Fragment_Start_Position, e.g., 12411</li>
 *         <li>Fragment_End_Position, e.g., 12460</li>
 *         <li>Fragment_Number, e.g., 3</li>
 *         <li>5'_Restriction_Site, e.g., DpnII</li>
 *         <li>3'_Restriction_Site, e.g., DpnII</li>
 *         <li>Length, e.g., 50</li>
 *         <li>5'_GC_Content, e.g., 0.620</li>
 *         <li>3'_GC_Content, e.g., 0.620</li>
 *         <li>5'_Repeat_Content, e.g., 0.000</li>
 *         <li>3'_Repeat_Content, e.g., 0.000 </li>
 *         <li>Selected, e.g., F</li>
 *         <li>5'_Probes, e.g., 0</li>
 *         <li>3'_Probes, e.g., 0</li>
 *     </ol>
 * </p>
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.2 (2020-11-18)
 */

public class DigestMap {
    private static final Logger logger = LoggerFactory.getLogger(DigestMap.class);
    /** Key: contig/chromosomes, with and without "chr". For instance, 'chr5' and '5' point to the
     * same {@link Chromosome2DigestArray} object.
     */
    private final Map<String, Chromosome2DigestArray> digestMap;

    public DigestMap(String digestFilePath) throws DiachromaticException {
        File f = new File(digestFilePath);
        if (! f.exists()) {
            throw new DiachromaticException(String.format("Could not find digest file at %s.", f.getAbsolutePath() ));
        }
        else {
            logger.trace("Found digest file at {}.",digestFilePath);
        }
        try ( BufferedReader br = new BufferedReader(new FileReader(digestFilePath))){
            Map<String, Chromosome2DigestArray> prelimMap = new HashMap<>();
            String line;
            while ((line=br.readLine())!=null) {
                if (line.startsWith("Chromosome")) continue; // the header line
                String[] fields = line.split("\t");
                if (fields.length!=Digest.TOTAL_NUMBER_OF_FIELDS) {
                    throw new DiachromaticException(String.format("Malformed line with %d fields (required: %d): %s",
                            fields.length,Digest.TOTAL_NUMBER_OF_FIELDS,line ));
                }
                Digest digest = new Digest(fields);
                String chromosome = digest.getChromosome();
                prelimMap.putIfAbsent(chromosome,new Chromosome2DigestArray());
                prelimMap.get(chromosome).addDigest(digest);
            }
            br.close();

            if (prelimMap.containsKey("M") && ! prelimMap.containsKey("MT")) {
                prelimMap.put("MT", prelimMap.get("M"));
            }
            if (prelimMap.containsKey("MT") && ! prelimMap.containsKey("M")) {
                prelimMap.put("M", prelimMap.get("MT"));
            }
            if (prelimMap.containsKey("chrM") && ! prelimMap.containsKey("chrMT")) {
                prelimMap.put("chrMT", prelimMap.get("chrM"));
            }
            if (prelimMap.containsKey("chrMT") && ! prelimMap.containsKey("chrM")) {
                prelimMap.put("chrM", prelimMap.get("chrMT"));
            }
            // In some cases, our data uses "chr5" and in others we see just "5".
            // The following adds some additional references to mitigate this issue
            for (Map.Entry<String, Chromosome2DigestArray> e : prelimMap.entrySet()) {
                if (e.getKey().startsWith("chr")) {
                    String newKey = e.getKey().substring(3);
                    prelimMap.put(newKey, e.getValue());
                } else {
                    String newKey = "chr" + e.getKey();
                    prelimMap.put(newKey, e.getValue());
                }
            }
            this.digestMap = Map.copyOf(prelimMap); // make immutable
        } catch (IOException e){
            throw new DiachromaticException(String.format("Could not parse %s: %s",digestFilePath,e.getMessage()));
        }
    }

    public Map<String, Chromosome2DigestArray> getDigestMap() {
        return digestMap;
    }


    /**
     * Get the pair of digests that correspond to the two positions defined by (chrom1,coord1) and (chrom2,coord2)
     * @param chrom1 chromosome of forward read
     * @param coord1 position on the chromosome of forward read
     * @param chrom2 chromosome of reverse read
     * @param coord2 position on the chromosome of reverse read
     * @return DigestPair object corresponding to these position.
     */
    public DigestPair getDigestPair(String chrom1, int coord1, String chrom2, int coord2) {
        // handle exception with unknown reference IDs see test
        Chromosome2DigestArray chrom2array1 = this.digestMap.get(chrom1);
        Digest digest1 = chrom2array1.getDigestAt(coord1);
        Chromosome2DigestArray chrom2array2 = this.digestMap.get(chrom2);
        Digest digest2 = chrom2array2.getDigestAt(coord2);
        return new DigestPair(digest1, digest2);
    }


    /**
     * This class stores all of the digests that are located on one chromosome (or scaffold). It additionally
     * stores an array of locations (the end positions) of each of the digests so that we can find the correct
     * digest given a position quicly using a binary search. The class is intended to be used with a map whose
     * key stores the name of the chromosome; the values of the map are objects of this class (one per chromosome).
     */
    static class Chromosome2DigestArray {
        /** List of the chromosomal positions of the digests on this chromosome. The end position is stored for each digest.*/
        private final ArrayList<Integer> coordArray;
        /** List of {@link Digest} objects corresponding to this chromosome. */
        private final ArrayList<Digest> digestArray;


        Chromosome2DigestArray() {
            coordArray = new ArrayList<>();
            digestArray = new ArrayList<>();
        }

        void addDigest(Digest digest) {
            digestArray.add(digest);
            coordArray.add(digest.getDigestEndPosition());
        }

        /**
         * Collections.binary search returns the index of a key for a list tyhat is sorted in ascending order.
         * If a key is not present, the method returns "(-(insertion point) - 1), whereby the insertion point is
         * defined as the point at which the key would be inserted into the list. The expression in the
         * "else" statement reverse this. We therefore get the right digest whether the search is made with
         * the actual END POSITION (which is stored in the array) or whether an INTERNAL POSITION is used.
         * @param coordinate
         * @return
         */
        public Digest getDigestAt(int coordinate){
            int index = Collections.binarySearch(this.coordArray, coordinate);
            Digest digest;
            if(0 <= index) {
                // coord1 is in the list and corresponds to digest end position
                digest = digestArray.get(index);
            } else {
                // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
                int i = (index+1)*(-1);
                //logger.trace("index: " + index + " i: " + i);
                digest = digestArray.get(i);
            }
            return digest;
        }

        public int getNumOfDigestsForChromosome() {
            return digestArray.size();
        }

    }

}

