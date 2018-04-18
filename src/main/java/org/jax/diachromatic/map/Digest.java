package org.jax.diachromatic.map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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

    private final int startpos;

    private final int endpos;

    private final int fragmentNumber;


    private String fivePrimeRestrictionSite;

    private String threePrimeRestrictionSite;

    private boolean active = false;


    public Digest(String[] fields) {
        chromosome=fields[0];
        startpos=Integer.parseInt(fields[1]);
        endpos=Integer.parseInt(fields[2]);
        fragmentNumber=Integer.parseInt(fields[3]);
        fivePrimeRestrictionSite=fields[4];
        threePrimeRestrictionSite=fields[5];
    }


    public String getChromosome() {
        return chromosome;
    }

    int getStartpos() {
        return startpos;
    }

    int getEndpos() {
        return endpos;
    }

    int getFragmentNumber() {
        return fragmentNumber;
    }

    public String getFivePrimeRestrictionSite() {
        return fivePrimeRestrictionSite;
    }

    public String getThreePrimeRestrictionSite() {
        return threePrimeRestrictionSite;
    }


    public int getSize() {
        return endpos - startpos + 1;
    }


    private void setActice() {
        this.active=true;
    }

    public boolean isActive() {
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
        startpos==other.startpos);
    }


    /**
     * Parse in the digest file (see {@link org.jax.diachromatic.digest.FragmentFactory} for details on file format).
     * The map has the chromosome as a key and a list of {@link Digest} objects on the chromosome as the value
     * @param digestFilePath path to the digest file created by {@link org.jax.diachromatic.command.DigestCommand}.
     * @return
     */
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
                    String fields[] = line.split("\t");
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
                String fields[] = line.split("\t");
                if (fields.length< 6) {
                    logger.fatal(String.format("Malformed line with %d fields (required: at least 6): %s",fields.length,line ));
                    System.exit(1); // todo throw exception
                }
                Digest dig = new Digest(fields);
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
                key += dig.getStartpos();
                key += "-";
                key += dig.getEndpos();
                if(activeDigestsFile != null && activeDigests.contains(key)) {
                    dig.setActice();
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
                startpos,
                endpos,
                fragmentNumber,
                fivePrimeRestrictionSite,
                threePrimeRestrictionSite);
    }
}
