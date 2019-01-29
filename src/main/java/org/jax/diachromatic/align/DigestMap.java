package org.jax.diachromatic.align;


import htsjdk.samtools.util.Log;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

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
 * @version 0.0.0 (2018-05-09)
 */

public class DigestMap {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    /** NEW VERSION */
    private final HashMap<String, Chromosome2DigestArray> digestMap;

    public DigestMap(String digestFilePath) throws DiachromaticException {
        this.digestMap = new HashMap<>();
        try {
            parseDigestFile(digestFilePath);
        } catch (IOException e){
            throw new DiachromaticException(String.format("Could not parse %s: %s",digestFilePath,e.getMessage()));
        }
    }



    public void parseDigestFile(String digestFilePath) throws IOException, DiachromaticException {
        File f = new File(digestFilePath);
        if (! f.exists()) {
            throw new DiachromaticException(String.format("Could not find digest file at %s", f.getAbsolutePath() ));
        }
        else {
            logger.trace("Found digest file at {}.",digestFilePath);
        }

        BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
        String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("Chromosome")) continue; // the header line
            String fields[] = line.split("\t");
            if (fields.length!=Digest.TOTAL_NUMBER_OF_FIELDS) {
                throw new DiachromaticException(String.format("Malformed line with %d fields (required: %d): %s",
                        Digest.TOTAL_NUMBER_OF_FIELDS,fields.length,line ));
            }
            Digest digest = new Digest(fields);
            String chromosome = digest.getChromosome();
            Integer digestEnd = digest.getDigestEndPosition();
            digestMap.putIfAbsent(chromosome,new Chromosome2DigestArray());
            digestMap.get(chromosome).addDigest(digest);
        }
        br.close();
        /*
        DO WE NEED THIS OR CAN/SHOULD WE ASSUME GOPHER WILL SORT?
        // sort position arrays for each chromosome and init activeState arrays
            // --------------------------------------------------------------------

            for (String key : digestMap.keySet()) {
                digestMap.get(key).finalizeArrayPair();
            }
         */
    }




    public DigestPair getDigestPair(String chrom1, int coord1, String chrom2, int coord2) {
        // handle exception with unknown reference IDs see test
        int index = Collections.binarySearch(this.digestMap.get(chrom1).coordArray, coord1);
        Digest digest1, digest2;
        if(0 <= index) {
            // coord1 is in the list and corresponds to digest end position
            digest1 = this.digestMap.get(chrom1).digestArray.get(index);
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            digest1 = this.digestMap.get(chrom1).digestArray.get(i);
        }
        int index2 = Collections.binarySearch(this.digestMap.get(chrom2).coordArray, coord2);
        if(0 <= index2) {
            // coord1 is in the list and corresponds to digest end position
            digest2 = this.digestMap.get(chrom2).digestArray.get(index2);
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index2+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            digest2 = this.digestMap.get(chrom2).digestArray.get(i);
        }
        return new DigestPair(digest1, digest2);
    }


    /**
     * This class stores all of the digests that are located on one chromosome (or scaffold). It additionally
     * stores an array of locations (the end positions) of each of the digests so that we can find the correct
     * digest given a position quicly using a binary search. The class is intended to be used with a map whose
     * key stores the name of the chromosome; the values of the map are objects of this class (one per chromosome).
     */
    private static class Chromosome2DigestArray {
        /** List of the chromosomal positions of the digests on this chromosome. The end position is stored for each digest.*/
        private ArrayList<Integer> coordArray;
        /** List of {@link Digest} objects corresponding to this chromosome. */
        private ArrayList<Digest> digestArray;
        /** Set of coordinates of all of the active digests. TODO do we need this???????????? */
        @Deprecated
        Set<Integer> activeStateCoordSet;

        Chromosome2DigestArray() {
            coordArray = new ArrayList<>();
            digestArray = new ArrayList<>();
            activeStateCoordSet = new HashSet<>();
        }

        void addDigest(Digest digest) {
            digestArray.add(digest);
            coordArray.add(digest.getDigestEndPosition());
            // TODO DO WE NEED THIS?????
//            if (digest.isSelected()) {
//                activeStateCoordSet.add(digest.getDigestEndPosition());
//            }
        }

        /*  REPLACED BY addDigest
        public void addCoord(Integer x) {
            coordArray.add(x);
        }*/

        public Integer getCoord(Integer index) {
            return coordArray.get(index);
        }

        /*  REPLACED BY addDigest
        public void addActiveStateCoord(Integer x) {
            activeStateCoordSet.add(x);
        }*/

        /* Koennen wir nicht annehmen dass die Digests sortiert sind? */
        @Deprecated
        public void finalizeArrayPair() {

            Collections.sort(coordArray);
// Note we now store the "activity" within the digest object.
//            for(int i = 0; i < coordArray.size(); i++) {
//                if(activeStateCoordSet.contains(coordArray.get(i))) {
//                    stateArray.add(1);
//                }
//                else {
//                    stateArray.add(0);
//                }
//            }
        }
    }

}

