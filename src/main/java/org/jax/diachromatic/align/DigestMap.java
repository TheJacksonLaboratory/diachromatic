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

    private HashMap<String,ArrayPair> digestMap = null; // key: chromosome name, value: pair of arrays (start coordinate of digest and active state)
    /** NEW VERSION */
    private HashMap<String,ArrayPair2> digestMap2;

    public DigestMap(String digestFilePath, String activeDigestsFile) throws DiachromaticException {

        try {

            // create hash map from the file for active digests
            Set<String> activeDigests = new HashSet<>();

            if (activeDigestsFile == null) {
                logger.trace(String.format("No file for active digests available. Will set all digests to inactive."));
            }
            else {
                logger.trace(String.format("File for active digests available. Reading file..."));

                File af = new File(activeDigestsFile);
                BufferedReader br = new BufferedReader(new FileReader(activeDigestsFile));
                String line;
                while ((line=br.readLine())!=null) {
                    String fields[] = line.split("\t");
                    if (fields.length < 3) {
                        logger.fatal(String.format("Malformed line with %d fields (required: at least 3): %s",fields.length,line ));
                        System.exit(1); // TODO: Add proper exception handling
                    }
                    String key=fields[0];
                    key += ":";
                    key += fields[2]; // end coordinate
                    if(activeDigestsFile != null) {
                        activeDigests.add(key);
                    }
                }
            }

            // read all digest coordinates and fill array digestEndPosition
            // ------------------------------------------------------------

            digestMap = new HashMap<String,ArrayPair>();

            File f = new File(digestFilePath);
            if (! f.exists()) {
                logger.error(String.format("Could not find digest file at %s", f.getAbsolutePath() ));
                System.exit(1);
            }
            else {
                logger.trace("File for all digests available. Reading file...");
            }

            BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
            String line;
            while ((line=br.readLine())!=null) {
                if (line.startsWith("Chromosome")) continue; // the header line
                String fields[] = line.split("\t");
                if (fields.length< 6) {
                    logger.fatal(String.format("Malformed line with %d fields (required: at least 6): %s",fields.length,line ));
                    System.exit(1); // todo throw exception
                }
                String chromosome = fields[0];
                Integer digestEnd = Integer.parseInt(fields[2]);
                if (!digestMap.containsKey(chromosome)) {
                    digestMap.put(chromosome,new ArrayPair());
                }

                digestMap.get(chromosome).addCoord(digestEnd);
                String key = chromosome;
                key += ":";
                key += digestEnd;
                if(activeDigests.contains(key))
                {
                    digestMap.get(chromosome).addActiveStateCoord(digestEnd);
                }

                /*
                 An additional file for active digests passed with -a option overwrites infos in columns 11 of the
                 digest file. Only if no file for active fragments is given, the information in column 11 is used.
                  */
                if(activeDigestsFile == null) {
                    if(fields[11].equals("T")) {
                        digestMap.get(chromosome).addActiveStateCoord(digestEnd);
                    }
                }
            }
            br.close();


            // sort position arrays for each chromosome and init activeState arrays
            // --------------------------------------------------------------------

            for (String key : digestMap.keySet()) {
                digestMap.get(key).finalizeArrayPair();
            }

        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1); // todo throw exception
        }
    }



    public void parseDigestFile(String digestFilePath) throws IOException, DiachromaticException {
        digestMap2 = new HashMap<>();

        File f = new File(digestFilePath);
        if (! f.exists()) {
            logger.error(String.format("Could not find digest file at %s", f.getAbsolutePath() ));
            System.exit(1);
        }
        else {
            logger.trace("File for all digests available. Reading file...");
        }

        BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
        String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("Chromosome")) continue; // the header line
            String fields[] = line.split("\t");
            if (fields.length!=Digest.TOTAL_NUMBER_OF_FIELDS) {
                logger.fatal(String.format("Malformed line with %d fields (required: %d): %s",
                        Digest.TOTAL_NUMBER_OF_FIELDS,fields.length,line ));
                System.exit(1); // todo throw exception
            }
            Digest digest = new Digest(fields);
            String chromosome = digest.getChromosome();
            Integer digestEnd = digest.getDigestEndPosition();
            digestMap2.putIfAbsent(chromosome,new ArrayPair2());
            digestMap2.get(chromosome).addDigest(digest);
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





    public DigestPair getDigestPair2(String chrom1, Integer coord1, String chrom2, Integer coord2) {

        // handle exception with unknown reference IDs see test
        int index = Collections.binarySearch(this.digestMap.get(chrom1).coordArray, coord1);
        String d1[] = new String[6];
        d1[0] = chrom1;
        Integer staCoord, endCoord;
        if(0 <= index) {
            // coord1 is in the list and corresponds to digest end position
            endCoord = this.digestMap.get(chrom1).coordArray.get(index);
            if(index == 0) {
                // this is the first digest
                staCoord = 1;
            } else {
                staCoord = this.digestMap.get(chrom1).coordArray.get(index-1) + 1;
            }
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            endCoord = this.digestMap.get(chrom1).coordArray.get(i);
            if(i == 0) {
                // this is the first digest
                staCoord = 0;
            } else {
                staCoord = this.digestMap.get(chrom1).coordArray.get(i-1) + 1;
            }
        }
        d1[1] = staCoord.toString();
        d1[2] = endCoord.toString();
        d1[3] = "42";
        d1[4] = "Dpn2";
        d1[5] = "Dpn2";
        /*
        Digest digest1 = new Digest(d1);
        if(this.digestMap.get(chrom1).activeStateCoordSet.contains(endCoord)) {
            digest1.setSelected();
        }
        index = Collections.binarySearch(this.digestMap.get(chrom2).coordArray, coord2);
        String d2[] = new String[6];
        d2[0] = chrom2;
        if(0 <= index) {
            // coord1 is in the list and corresponds to digest end position
            endCoord = this.digestMap.get(chrom2).coordArray.get(index);
            if(index == 0) {
                // this is the first digest
                staCoord = 1;
            } else {
                staCoord = this.digestMap.get(chrom2).coordArray.get(index-1) + 1;
            }
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            endCoord = this.digestMap.get(chrom2).coordArray.get(i);
            if(i == 0) {
                // this is the first digest
                staCoord = 0;
            } else {
                staCoord = this.digestMap.get(chrom2).coordArray.get(i-1) + 1;
            }
        }
        d2[1] = staCoord.toString();
        d2[2] = endCoord.toString();
        d2[3] = "43";
        d2[4] = "Dpn2";
        d2[5] = "Dpn2";
        Digest digest2 = new Digest(d2);
        if(this.digestMap.get(chrom2).activeStateCoordSet.contains(endCoord)) {
            digest2.setSelected();
        }
        return new DigestPair(digest1, digest2);
        */
        return null;
    }


    public DigestPair getDigestPairNewAndImproved(String chrom1, int coord1, String chrom2, int coord2) {
        // handle exception with unknown reference IDs see test
        int index = Collections.binarySearch(this.digestMap2.get(chrom1).coordArray, coord1);
        Digest digest1, digest2;
        if(0 <= index) {
            // coord1 is in the list and corresponds to digest end position
            digest1 = this.digestMap2.get(chrom1).digestArray.get(index);
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            digest1 = this.digestMap2.get(chrom1).digestArray.get(i);
        }
        int index2 = Collections.binarySearch(this.digestMap2.get(chrom2).coordArray, coord2);
        if(0 <= index2) {
            // coord1 is in the list and corresponds to digest end position
            digest2 = this.digestMap2.get(chrom2).digestArray.get(index2);
        } else {
            // coord1 is not in the list and would be inserted at i=(index+1)*(-1)
            int i = (index2+1)*(-1);
            //logger.trace("index: " + index + " i: " + i);
            digest2 = this.digestMap2.get(chrom2).digestArray.get(i);
        }
        return new DigestPair(digest1, digest2);
    }



    private class ArrayPair {

        private ArrayList<Integer> coordArray;
        private ArrayList<Integer> stateArray;
        Set<Integer> activeStateCoordSet;

        ArrayPair() {
            coordArray = new ArrayList<Integer>();
            stateArray = new ArrayList<Integer>();
            activeStateCoordSet = new HashSet<>();
        }

        public void addCoord(Integer x) {
            coordArray.add(x);
        }

        public Integer getCoord(Integer index) {
            return coordArray.get(index);
        }

        public void addActiveStateCoord(Integer x) {
            activeStateCoordSet.add(x);
        }

        public void finalizeArrayPair() {

            Collections.sort(coordArray);

            for(int i = 0; i < coordArray.size(); i++) {
                if(activeStateCoordSet.contains(coordArray.get(i))) {
                    stateArray.add(1);
                }
                else {
                    stateArray.add(0);
                }
            }
        }
    }



    private static class ArrayPair2 {

        private ArrayList<Integer> coordArray;
        private ArrayList<Integer> stateArray;
        private ArrayList<Digest> digestArray;
        Set<Integer> activeStateCoordSet;

        ArrayPair2() {
            coordArray = new ArrayList<>();
            stateArray = new ArrayList<>();
            digestArray = new ArrayList<>();
            activeStateCoordSet = new HashSet<>();
        }

        void addDigest(Digest digest) {
            digestArray.add(digest);
            coordArray.add(digest.getDigestEndPosition());
            if (digest.isSelected()) {
                activeStateCoordSet.add(digest.getDigestEndPosition());
            }
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

            for(int i = 0; i < coordArray.size(); i++) {
                if(activeStateCoordSet.contains(coordArray.get(i))) {
                    stateArray.add(1);
                }
                else {
                    stateArray.add(0);
                }
            }
        }
    }

}

