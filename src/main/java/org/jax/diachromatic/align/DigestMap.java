package org.jax.diachromatic.align;


import htsjdk.samtools.SAMRecord;
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
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.0 (2018-05-09)
 */

public class DigestMap {
    private static final Logger logger = LogManager.getLogger();
    private static final htsjdk.samtools.util.Log log = Log.getInstance(Aligner.class);

    private HashMap<String,ArrayPair> digestMap = null; // key: chromosome name, value: pair of arrays (start coordinate of digest and active state)

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
                    key += fields[1];
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
                Integer digestStart = Integer.parseInt(fields[1]);
                if (!digestMap.containsKey(chromosome)) {
                    digestMap.put(chromosome,new ArrayPair());
                }
                digestMap.get(chromosome).addCoord(digestStart);
                String key = chromosome;
                key += ":";
                key += digestStart;
                if(activeDigests.contains(key))
                {
                    digestMap.get(chromosome).addActiveStateCoord(digestStart);
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

    public DigestPair getDigestPair2(String chrom1, Integer coord1,String chrom2, Integer coord2) {
        int index = Collections.binarySearch(this.digestMap.get(chrom1).coordArray, coord1);
        String d1[] = new String[6];
        d1[0] = chrom1;
        d1[1] = this.digestMap.get(chrom1).coordArray.get(index).toString();
        d1[2] = this.digestMap.get(chrom1).coordArray.get(index+1).toString();
        d1[3] = "42";
        d1[4] = "Dpn2";
        d1[5] = "Dpn2";
        Digest digest1 = new Digest(d1);
        index = Collections.binarySearch(this.digestMap.get(chrom2).coordArray, coord2);
        String d2[] = new String[6];
        d2[0] = chrom1;
        d2[1] = this.digestMap.get(chrom2).coordArray.get(index).toString();
        d2[2] = this.digestMap.get(chrom2).coordArray.get(index+1).toString();
        d2[3] = "42";
        d2[4] = "Dpn2";
        d2[5] = "Dpn2";
        Digest digest2 = new Digest(d2);
        return new DigestPair(digest1, digest2);
    }



    private class ArrayPair {

        private ArrayList coordArray;
        private ArrayList stateArray;
        Set<Integer> activeStateCoordSet;

        ArrayPair() {
            coordArray = new ArrayList();
            stateArray = new ArrayList();
            activeStateCoordSet = new HashSet<>();
        }

        public void addCoord(Integer x) {
            coordArray.add(x);
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

}

