package org.jax.diachromatic.digest;

import com.sun.org.apache.regexp.internal.RE;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @version 0.2.3 (2017-11-25)
 */
public class RestrictionEnzyme implements Serializable {
    private static final Logger logger = LogManager.getLogger();
    /** serialization version ID */
    static final long serialVersionUID = 2L;
    /** A name of a restriction enzyme, something like HindIII */
    private String name;
    /** A representation of the cutting site of the enzyme, whereby "^" stands for cut here.
     * For instance, A^AGCTT is the cutting site for HindIII.
     */
    private String site;

    /** Same as site but without the caret symbol */
    private String plainSite;
    /** The offset of the cutting site in this restriction enzyme. For instancen the offset for ^GATC is 0 and the
     * offset for A^AGCTT is 1.
     */
    private Integer offset=null;

    public RestrictionEnzyme(String n, String s) {
        name=n;
        site=s;
        this.offset=site.indexOf('^');
        if (offset<0) {
            logger.error(String.format("Malformed site pattern for enyze %s (%s)",name,site)); /* Should never happen!*/
        }
        plainSite=site;
        int i= site.indexOf('^');
        if (i>=0) {
            plainSite=site.substring(0,i)+site.substring(i+1);
        }
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getSite() {
        return site;
    }

    public String getPlainSite() { return  plainSite; }

    public void setSite(String site) {
        this.site = site;
    }

    public String getLabel() {
        return String.format("%s: %s",getName(),getSite());
    }

    public Integer getOffset() { return this.offset; }


    public static List<RestrictionEnzyme> parseRestrictionEnzymesFromFile(String path) {
        List<RestrictionEnzyme> reList=new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(path));
            String line=null;
            while ((line=br.readLine())!=null) {
                System.out.println();
                if (line.startsWith("#"))
                    continue; // comment
                String A[] = line.split("\\s+");
                if (A.length!=2) {
                    logger.error(String.format("Malformed restrtiction enzyme line-should have 2 fields but had %d (%s)",A.length,line));
                    continue;
                }
                RestrictionEnzyme re = new RestrictionEnzyme(A[0].trim(),A[1].trim());
                reList.add(re);

            }
            br.close();
        } catch (IOException e) {
            logger.fatal(String.format("Could not read restriction enzymes from %s",path));
        }
        return reList;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj==null) return false;
        if (! (obj instanceof RestrictionEnzyme)) return false;
        RestrictionEnzyme other = (RestrictionEnzyme) obj;
        return (name.equals(other.name) && site.equals(other.site));
    }


}
