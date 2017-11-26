package org.jax.diachromatic.truncation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.IOException;

public class FastQRecord {
    private static final Logger logger = LogManager.getLogger();
    private final String name;
    private final String name2;
    private  String sequence=null;
    private  String quality=null;

    private static String ligationSequence=null;
    private static String restrictionSequence=null;


    public static void setLigationSequence(String seq) {
        ligationSequence=seq;
    }

    public static void setRestrictionSequence(String seq) {
        restrictionSequence=seq;
    }


    public FastQRecord(String [] lines) {
        assert lines.length==4; // should never happen, by design we only pass 4 entry arrays, TODO remove if all is ok
        name=lines[0];
        sequence=lines[1];
        // note lines[2] should be just "."
        name2=lines[2];
        quality=lines[3];
        if (! name.startsWith("@")) {
            logger.fatal(String.format("Malformed FASTQ name line (did not start with @ symbol): %s",name ));
            System.exit(1); // to do throw exception
        }
    }

    /**
     * truncates sequence and quality lines if the {@link #ligationSequence} is found, and returns true, otherwise
     * returns false
     * @return true if truncation was performed.
     */
    public boolean truncateIfLigationSiteFound() {
        int i = sequence.indexOf(FastQRecord.ligationSequence);
        if (i<0) {
            return false;  // we did not find the ligation sequence
        }
        sequence=sequence.substring(0,i) + FastQRecord.restrictionSequence;
        int len=sequence.length();
        quality=quality.substring(0,len);
        return true;
    }


    public void writeToStream(BufferedWriter out) throws IOException {
        out.write(name+"\n");
        out.write(sequence + "\n");
        out.write(name2 +"\n");
        out.write(quality+"\n");
    }


}
