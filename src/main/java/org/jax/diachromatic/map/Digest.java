package org.jax.diachromatic.map;

import com.google.common.collect.ImmutableList;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

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


    /**
     * Parse in the digest file (see {@link org.jax.diachromatic.digest.FragmentFactory} for details on file format).
     * @param digestFilePath
     * @return
     */
    public static ImmutableList<Digest> readDigests(String digestFilePath) {
        ImmutableList.Builder<Digest>  builder = new ImmutableList.Builder();
        try {
            BufferedReader br = new BufferedReader(new FileReader(digestFilePath));
            String line=null;
            while ((line=br.readLine())!=null) {
                System.out.println(line);
                if (line.startsWith("Chromosome")) continue; // the header line
                String fields[] = line.split("\t");
                if (fields.length!= 6) {
                    logger.fatal(String.format("Malformed line with %d fields (required: 6): %s",fields.length,line ));
                    System.exit(1); // todo throw exception
                }
                Digest dig = new Digest(fields);
                builder.add(dig);
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1); // todo throw exception
        }
        return builder.build();
    }
}
