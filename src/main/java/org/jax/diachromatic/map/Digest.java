package org.jax.diachromatic.map;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Digest {
    private static final Logger logger = LogManager.getLogger();
    private final String chromosome;

    private final int startpos;

    private final int endpos;

    private final int fragmentNumber;


    private String fivePrimeRestrictionSite;

    private String threePrimeRestrictionSite;


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

    public int getStartpos() {
        return startpos;
    }

    public int getEndpos() {
        return endpos;
    }

    public int getFragmentNumber() {
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


    /** Note -- by the way these objects are created in this program, it is sufficient to check whether
     * the chromosome and the start position are equal in order to know whether the objects are equal.
     * @param o
     * @return
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
     * @param digestFilePath
     * @return
     */
    public static Map<String,List<Digest>> readDigests(String digestFilePath) {
        Map<String,List<Digest>> map = new HashMap<>();
        try {
            File f = new File(digestFilePath);
            if (! f.exists()) {
                logger.error(String.format("Could not find digest file at ", f.getAbsolutePath() ));
                System.exit(1);
            }

            BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
            String line=null;
            while ((line=br.readLine())!=null) {
                //System.out.println(line);
                if (line.startsWith("Chromosome")) continue; // the header line
                String fields[] = line.split("\t");
                if (fields.length!= 6) {
                    logger.fatal(String.format("Malformed line with %d fields (required: 6): %s",fields.length,line ));
                    System.exit(1); // todo throw exception
                }
                Digest dig = new Digest(fields);
                String chrom = dig.getChromosome();
                List<Digest> dlist=null;
                if (map.containsKey(chrom)) {
                    dlist=map.get(chrom);
                } else {
                    dlist=new ArrayList<>();
                    map.put(chrom,dlist);
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
}
