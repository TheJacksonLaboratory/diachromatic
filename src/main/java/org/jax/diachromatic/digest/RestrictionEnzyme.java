package org.jax.diachromatic.digest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Note that this class was taken originally from VPV and extended here to read files. TODO
 * Consider putting this version of the class back into VPV.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @version 0.2.4 (2017-11-25)
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

    /**
     * Dangling end sequence, i.e. the ends of un-ligated sequences. This is the sequence between the cutting site and
     * the last nucleotide of the recognition motif, e.g.
     *
     * for HindIII the recognition motif is A^AGCTT and the dangling end sequence is AGCTT
     *
     * or for DpnII the recognition motif is ^GATC and the dangling end sequence is GATC.
     *
     * The dangling end sequence is independent of the sticky ends option,
     * i.e. whether a fill-in of the sticky ends was performed for Hi-C or not.
     */
    private String danglingEndSeq;

    /** The offset of the cutting site in this restriction enzyme. For instancen the offset for ^GATC is 0 and the
     * offset for A^AGCTT is 1.
     */
    private Integer offset=null;

    public RestrictionEnzyme(String n, String s) {
        name = n;
        site = s;
        this.offset = site.indexOf('^');
        if (offset < 0) {
            // Should never happen!
            logger.error(String.format("Malformed site pattern for enyze %s (%s)", name, site));
        }
        plainSite = site;
        int i = site.indexOf('^');
        if (i >= 0) {
            plainSite = site.substring(0, i) + site.substring(i + 1);
        }
        int center = plainSite.length() / 2;
        if (offset <= center) {
            this.danglingEndSeq = plainSite.substring(offset, plainSite.length());
        } else {
            this.danglingEndSeq = plainSite.substring(plainSite.length() - offset, plainSite.length());
        }
    }

    public String getName() {
        return name;
    }

    public String getSite() {
        return site;
    }

    public String getPlainSite() { return  plainSite; }

    public String getLabel() {
        return String.format("%s: %s",getName(),getSite());
    }

    public Integer getOffset() { return this.offset; }

    public  String getDanglingEndSequence() {
        return this.danglingEndSeq;
    }


    public static List<RestrictionEnzyme> parseRestrictionEnzymes() {
        List<RestrictionEnzyme> reList = new ArrayList<>();
        ClassLoader classLoader = RestrictionEnzyme.class.getClassLoader();
        InputStream is = classLoader.getResourceAsStream("data/enzymelist.tab");
        try {
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String line;
            while ((line=br.readLine())!=null) {
                if (line.startsWith("#"))
                    continue; // comment
                String[] A = line.split("\\s+");
                if (A.length!=2) {
                    logger.error(String.format("Malformed restrtiction enzyme line-should have 2 fields but had %d (%s)",A.length,line));
                    continue;
                }
                RestrictionEnzyme re = new RestrictionEnzyme(A[0].trim(),A[1].trim());
                reList.add(re);
            }
            br.close();
        } catch (IOException e) {
            logger.fatal("Could not read restriction enzymes from stream");
            logger.fatal(e);
            System.exit(1);
        }
        logger.trace(String.format("Got %d enzyme definitions.", reList.size()));
        return reList;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) return false;
        if (!(obj instanceof RestrictionEnzyme)) return false;
        RestrictionEnzyme other = (RestrictionEnzyme) obj;
        return (name.equals(other.name) && site.equals(other.site));
    }


}
