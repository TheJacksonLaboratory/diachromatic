package org.jax.diachromatic.map;

import htsjdk.samtools.SAMRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.util.Pair;

import java.util.*;

/**
 * This class is used to capture a valid ditag and to check for duplicates. A ditag is defined by
 * <ol>
 *     <li>Forward read starting position (sonication cut site)</li>
 * <li>Forward read orientation</li>
 * <li>Reverse read starting position (sonication cut site)</li>
 * <li>Reverse read orientation</li>
 * </ol>
 * Note that to termine the start position,  SAM format reports the start position of the sequence as it appears in the
 * reference genome without regarding the strand of the read.(not the first base encountered
 * in the sequencing read). Therefore, we determine the outermost 5' and 3' ends on the reads.
 * Determine the strand to which the read mapped by analysing the SAM bitwise flag
  * and thus whether the sonication cut site maps to the start or end of the
  * genomic reference sequence
 */
public class DiTag {
    private static final Logger logger = LogManager.getLogger();
    private Character chrom1;
    private byte strand1;
    private Character chrom2;
    private byte strand2;
    private int sonicationStartF;
    private int sonicationStartR;

    private static Character currentChar=0;

    private static Map<Character,String> char2stringChrom=new HashMap<>();
    private static Map<String,Character> string2char=new HashMap<>();

    private static Set<DiTag> ditagset=new HashSet();

    private DiTag (){
        // disallow contructor
    }

    private DiTag(Tag t1, Tag t2) {
        chrom1=t1.chrom;
        strand1=t1.strand;
        sonicationStartF=t1.sonicationStart;
        chrom2=t2.chrom;
        strand2=t2.strand;
        sonicationStartR=t2.sonicationStart;
    }


    /** We represent the chromosomes as (two-byte) Character objects to save space). */
    private static Character getChromosomeCharacter(String c) {
        Character chromchar;
        if ( string2char.containsKey(c)   ) {
            chromchar=string2char.get(c);
        } else {
            chromchar=currentChar;
            currentChar++;
            string2char.put(c,chromchar);
            char2stringChrom.put(chromchar,c);
        }
        return chromchar;
    }


    /** Note that which read of a pair is sequenced first is by chance -- we therefore first
     * record the relevant data for each read as a Tag, sort the Tags, and then create a {@link DiTag} object
     * from which we will determine if a given read pais is a duplicate or not.
     */
    static class Tag {
         Character chrom;
         byte strand;
        int sonicationStart;
    }


    /**
     * See explanation for {@link Tag}. We compare the order of two {@link Tag} objects with this function
     * prior to combining a pair of {@link Tag}'s into a {@link DiTag}.
     * @param t1
     * @param t2
     * @return
     */
    public static int compare(Tag t1, Tag t2)
    {
        if (! t1.chrom.equals(t2.chrom)) return t1.chrom.compareTo(t2.chrom);
        if (! (t1.strand == t2.strand) ) return t1.strand - t2.strand;
        return t1.sonicationStart - t2.sonicationStart;
    }

    public static boolean isDuplicate(ReadPair readpair) {

        Tag tag1 = new Tag();
        Tag tag2 = new Tag();
        tag1.chrom =getChromosomeCharacter(readpair.forward().getReferenceName());
        tag2.chrom =getChromosomeCharacter(readpair.reverse().getReferenceName());
        if (readpair.forward().getReadNegativeStrandFlag()) {
            tag1.strand ='-';
            tag1.sonicationStart =readpair.forward().getAlignmentEnd();
        } else {
            tag1.strand ='+';
            tag1.sonicationStart =readpair.forward().getAlignmentStart();
        }
        if (readpair.reverse().getReadNegativeStrandFlag()) {
            tag2.strand ='+';
            tag2.sonicationStart =readpair.reverse().getAlignmentEnd();
        } else {
            tag2.strand ='+';
            tag2.sonicationStart =readpair.reverse().getAlignmentStart();;
        }
        DiTag ditag=null;
        if (compare(tag1,tag2)>=0) {
             ditag = new DiTag(tag1,tag2);
        } else {
             ditag = new DiTag(tag2,tag1);
        }
        logger.trace(String.format("ditag %s",ditag.chrom1));
        if (ditagset.contains(ditag))  {
            return true; // we have already seen this readpair, it is a duplicate.
        } else {
            ditagset.add(ditag); // never before seen, add to HashSet
            return false;
        }
    }
}


