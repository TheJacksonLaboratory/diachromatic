package org.jax.diachromatic.map;

import htsjdk.samtools.SAMRecord;
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

    private DiTag(Tag[] ar) {
        chrom1=ar[0].chrom;
        strand1=ar[0].strand;
        sonicationStartF=ar[0].sonicationStart;
        chrom2=ar[1].chrom;
        strand2=ar[1].strand;
        sonicationStartR=ar[1].sonicationStart;
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



    static class Tag {
         Character chrom;
         byte strand;
        int sonicationStart;
    }

    static class TagComparator implements Comparator<Tag>
    {
        public int compare(Tag t1, Tag t2)
        {
            if (! t1.chrom.equals(t2.chrom)) return t1.chrom.compareTo(t2.chrom);
            if (! (t1.strand == t2.strand) ) return t1.strand - t2.strand;
            return t1.sonicationStart - t2.sonicationStart;
        }
    }



    public static boolean isDuplicate(Pair<SAMRecord,SAMRecord> readpair) {

        Tag tag1 = new Tag();
        Tag tag2 = new Tag();
        tag1.chrom =getChromosomeCharacter(readpair.first.getReferenceName());
        tag2.chrom =getChromosomeCharacter(readpair.second.getReferenceName());
        if (readpair.first.getReadNegativeStrandFlag()) {
//            ditag.strand='-';
//            ditag.sonicationStartF=readpair.first.getAlignmentEnd();
            tag1.strand ='-';
            tag1.sonicationStart =readpair.first.getAlignmentEnd();

        } else {
//            ditag.strand='+';
//            ditag.sonicationStartF=readpair.first.getAlignmentStart();
            tag1.strand ='+';
            tag1.sonicationStart =readpair.first.getAlignmentStart();
        }
        if (readpair.second.getReadNegativeStrandFlag()) {
//            ditag.strand2='-';
//            ditag.sonicationStartR=readpair.second.getAlignmentEnd();
            tag2.strand ='+';
            tag2.sonicationStart =readpair.second.getAlignmentEnd();

        } else {
//            ditag.strand2='+';
//            ditag.sonicationStartR=readpair.second.getAlignmentStart();
            tag2.strand ='+';
            tag2.sonicationStart =readpair.second.getAlignmentStart();;
        }
        Tag[] ar={tag1,tag2};
        Arrays.sort(ar,new TagComparator());

        DiTag ditag = new DiTag(ar);
        if (ditagset.contains(ditag)) return true;
        ditagset.add(ditag);
        return false;
//        ditag.chrom=getChromosomeCharacter(readpair.first.getReferenceName());
//        ditag.chrom2=getChromosomeCharacter(readpair.second.getReferenceName());

    }
}


/*
#############################
#Subroutine 'ditag_labeller':
#input is a paired read from a SAM file and returns a di-tag label comprising the
#chromsome name, start position and strand of each read in the pair
#A di-tag is defined by:
#Forward read starting position (sonication cut site)
#Forward read orientation
#Reverse read starting position (sonication cut site)
#Reverse read orientation
sub ditag_labeller {

    my ( $readF, $readR ) = @_;

    my $readF_strand = ( split( /\t/, $readF ) )[1];
    my $readF_csome  = ( split( /\t/, $readF ) )[2];
    my $readF_start  = ( split( /\t/, $readF ) )[3];    #Most upstream base on reference genome
    my $readF_seq    = ( split( /\t/, $readF ) )[9];
    my $readF_sonic;

    my $readR_strand = ( split( /\t/, $readR ) )[1];
    my $readR_csome  = ( split( /\t/, $readR ) )[2];
    my $readR_start  = ( split( /\t/, $readR ) )[3];    #Most upstream base on reference genome
    my $readR_seq    = ( split( /\t/, $readR ) )[9];
    my $readR_sonic;

    #
    if ( $readF_strand & 0x10 ) {
        $readF_strand = '-';
        $readF_sonic  = $readF_start + length($readF_seq) - 1;
    } else {
        $readF_strand = '+';
        $readF_sonic  = $readF_start;
    }

    if ( $readR_strand & 0x10 ) {
        $readR_strand = '-';
        $readR_sonic  = $readR_start + length($readR_seq) - 1;
    } else {
        $readR_strand = '+';
        $readR_sonic  = $readR_start;
    }

    my $labelF = "$readF_csome\t$readF_sonic\t$readF_strand";
    my $labelR = "$readR_csome\t$readR_sonic\t$readR_strand";
    my $ditag_label;

    #Which read of the read pair is sequenced first is random (i.e. expected to be 50:50),
    #and is not of biological significance.
    if ( ( $labelF cmp $labelR ) == 1 ) {
        $ditag_label = "$labelR\t$labelF";
    } else {
        $ditag_label = "$labelF\t$labelR";
    }
}

 */
