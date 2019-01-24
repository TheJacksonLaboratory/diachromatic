package org.jax.diachromatic.align;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * A class to represent a restriction fragment resulting from an in silico digestion of the genome.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.3 (2018-01-07)
 */
public class Digest {
    private static final Logger logger = LogManager.getLogger();
    private final String chromosome;

    private final int digestStartPosition;

    private final int digestEndPosition;
    /** We keep track of the digests by numbering them 1 to N, and use the {@link DigestMap} class to store and find them
     * according to this number.*/
    private final int digesttNumber;
    /** Length of the Digest (equal to {@link #digestEndPosition} = {@link #digestStartPosition} + 1). */
    private final int digestLength;


    private String fivePrimeRestrictionSite;

    private String threePrimeRestrictionSite;

    private final double five_prime_GC;

    private final double three_prime_GC;

    private final double five_prime_repeat;

    private final double three_prime_repeat;

    private final double five_prime_probe_count;

    private final double three_prime_probe_count;


    /** If true, then this digest has been selected for enrichment by a capture probe. */
    private boolean active = false;


    private final static int CHROMOSOME_INDEX=0;
    private final static int DIGEST_START_POSITION_INDEX=1;
    private final static int DIGEST_END_POSITION_INDEX=2;
    private final static int DIGEST_NUMBER_INDEX=3;
    private final static int FIVE_PRIME_RESTRICTION_SITE_INDEX=4;
    private final static int THREE_PRIME_RESTRICTION_SITE_INDEX=5;
    private final static int DIGEST_LENGTH_INDEX =6;
    private final static int FIVE_PRIME_GC_CONTENT_INDEX=7;
    private final static int THREE_PRIME_GC_CONTENT_INDEX=8;
    private final static int FIVE_PRIME_REPEAT_CONTENT_INDEX=9;
    private final static int THREE_PRIME_REPEAT_CONTENT_INDEX=10;
    private final static int SELECTED_INDEX=11;
    private final static int FIVE_PRIME_PROBE_COUNT_INDEX=12;
    private final static int THREE_PRIME_PROBE_COUNT_INDEX=13;
    /** total number of fields in the GOPHER digest file (separated by tabs). */
    public final static int TOTAL_NUMBER_OF_FIELDS=14;


    public Digest(String[] fields) throws DiachromaticException{
        if (fields.length != TOTAL_NUMBER_OF_FIELDS) {
            throw new DiachromaticException(String.format("Incorrect number of fields in digest file line: %d (%s)",
                    fields.length,
                    String.join(";",fields)));
        }
        chromosome=fields[CHROMOSOME_INDEX];
        digestStartPosition =Integer.parseInt(fields[DIGEST_START_POSITION_INDEX]);
        digestEndPosition =Integer.parseInt(fields[DIGEST_END_POSITION_INDEX]);
        digesttNumber =Integer.parseInt(fields[DIGEST_NUMBER_INDEX]);
        fivePrimeRestrictionSite=fields[FIVE_PRIME_RESTRICTION_SITE_INDEX];
        threePrimeRestrictionSite=fields[THREE_PRIME_RESTRICTION_SITE_INDEX];
        digestLength=Integer.parseInt(fields[DIGEST_LENGTH_INDEX]);
        five_prime_GC=Double.parseDouble(fields[FIVE_PRIME_GC_CONTENT_INDEX]);
        three_prime_GC=Double.parseDouble(fields[THREE_PRIME_GC_CONTENT_INDEX]);
        five_prime_repeat=Double.parseDouble(fields[FIVE_PRIME_REPEAT_CONTENT_INDEX]);
        three_prime_repeat=Double.parseDouble(fields[THREE_PRIME_REPEAT_CONTENT_INDEX]);
        if (fields[SELECTED_INDEX].equals("F")) {
            active=false;
        } else if (fields[SELECTED_INDEX].equals("T")) {
            active=true;
        } else {
            throw new DiachromaticException(String.format("Malformed selected field (%s) digest file line: (%s)",
                    fields[SELECTED_INDEX],
                    String.join(";",fields)));
        }
        five_prime_probe_count =Integer.parseInt(fields[FIVE_PRIME_PROBE_COUNT_INDEX]);
        three_prime_probe_count =Integer.parseInt(fields[THREE_PRIME_PROBE_COUNT_INDEX]);
    }







    public String getChromosome() {
        return chromosome;
    }

    int getDigestStartPosition() {
        return digestStartPosition;
    }

    int getDigestEndPosition() {
        return digestEndPosition;
    }

    int getDigesttNumber() {
        return digesttNumber;
    }

    public String getFivePrimeRestrictionSite() {
        return fivePrimeRestrictionSite;
    }

    public String getThreePrimeRestrictionSite() {
        return threePrimeRestrictionSite;
    }


    public int getSize() {
        return digestEndPosition - digestStartPosition + 1;
    }


    public void setSelected() {
        this.active=true;
    }

    public boolean isSelected() {
        return this.active;
    }

    /** Note -- by the way these objects are created in this program, it is sufficient to check whether
     * the chromosome and the start position are equal in order to know whether the objects are equal.
     * @param o the Object being compared with this.
     * @return true if o and this are equal
     */
    @Override
    public boolean equals(Object o) {
        if (o==null) return false;
        if (! (o instanceof Digest) ) return false;
        Digest other = (Digest) o;
        return (chromosome.equals(other.chromosome) &&
        digestStartPosition ==other.digestStartPosition);
    }


    /**
     * TODO do we need this??
     * Parse in the digest file for details on file format).
     * The align has the chromosome as a key and a list of {@link Digest} objects on the chromosome as the value
     * @param digestFilePath path to the digest file.
     * @return
     */
    @Deprecated
    public static Map<String,List<Digest>> readDigests(String digestFilePath, String activeDigestsFile) {
        Map<String,List<Digest>> map = new HashMap<>();
        try {

            Set<String> activeDigests = new HashSet<>();

            File f = new File(digestFilePath);
            if (! f.exists()) {
                logger.error(String.format("Could not find digest file at %s", f.getAbsolutePath() ));
                System.exit(1);
            }

            if (activeDigestsFile == null) {
                logger.trace(String.format("No file for active digests available. Will set all digests to inactive."));
            }
            else {
                logger.trace(String.format("File for active digests available. Reading file..."));

                File af = new File(activeDigestsFile);
                BufferedReader br = new BufferedReader(new FileReader(activeDigestsFile));
                String line;
                while ((line=br.readLine())!=null) {
                    String[] fields = line.split("\t");
                    if (fields.length < 3) {
                        logger.fatal(String.format("Malformed line with %d fields (required: at least 3): %s",fields.length,line ));
                        System.exit(1); // TODO: Add proper exception handling
                    }
                    String key=fields[0];
                    key += ":";
                    key += fields[1];
                    key += "-";
                    key += fields[2];
                    if(activeDigestsFile != null) {
                        activeDigests.add(key);
                    }
                }
            }

            BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
            String line;
            while ((line=br.readLine())!=null) {
                //System.out.println(line);
                if (line.startsWith("Chromosome")) continue; // the header line
                String[] fields = line.split("\t");
                if (fields.length< 6) {
                    logger.fatal(String.format("Malformed line with %d fields (required: at least 6): %s",fields.length,line ));
                    System.exit(1); // todo throw exception
                }
                Digest dig;
                try {
                    dig = new Digest(fields);
                } catch (DiachromaticException e) {
                    e.printStackTrace();
                    continue;
                }
                String chrom = dig.getChromosome();

                List<Digest> dlist;
                if (map.containsKey(chrom)) {
                    dlist=map.get(chrom);
                } else {
                    dlist=new ArrayList<>();
                    map.put(chrom,dlist);
                }

                // check if digest is active
                String key = dig.getChromosome();
                key += ":";
                key += dig.getDigestStartPosition();
                key += "-";
                key += dig.getDigestEndPosition();
                if(activeDigestsFile != null && activeDigests.contains(key)) {
                    dig.setSelected();
                }
                dlist.add(dig);
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1); // todo throw exception
        }
        return map;
    }

    @Override
    public String toString() {
        return String.format("Digest at %s:%d-%d [frag. %d;%s/%s]",
                chromosome,
                digestStartPosition,
                digestEndPosition,
                digesttNumber,
                fivePrimeRestrictionSite,
                threePrimeRestrictionSite);
    }
}
