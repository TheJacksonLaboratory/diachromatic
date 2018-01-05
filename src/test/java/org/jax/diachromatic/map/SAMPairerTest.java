package org.jax.diachromatic.map;

import htsjdk.samtools.SAMRecord;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;
import org.jax.diachromatic.util.Pair;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

public class SAMPairerTest {

    private static int digestcounter=0;
    private static String restrictionsite=null;
    private static Map<String,List<Digest>> digestmap;
    private static SAMPairer sampairer;
    private static String sam1;
    private static String sam2;

    @BeforeClass
    public  static void init() {
        ClassLoader classLoader = SAMPairerTest.class.getClassLoader();
        sam1 = classLoader.getResource("data/sam/forwardtest.sam").getFile();
        sam2 = classLoader.getResource("data/sam/reversetest.sam").getFile();
        // make fake digests for the reads we will test
        restrictionsite="HindIII";
        digestmap=new HashMap<>();
        List<Digest> digests = makeFakeDigestList("chr16",31526917);
        digestmap.put("chr16",digests);
        digests=makeFakeDigestList("chr11",92316468);
        digestmap.put("chr11",digests);

    }


    /**
     * We require a digest list for the {@link SAMPairer} constructor. We make a "fake" digest list for
     * simplicity that will allow us to perform testing of various functionalities.
     * @param chrom
     * @param pos
     * @return
     */
    private static List<Digest> makeFakeDigestList(String chrom,Integer ...pos) {
        List<Digest> dlist = new ArrayList<>();
        for (Integer p : pos) {
            Digest d = makeFakeDigest(chrom,p);
            dlist.add(d);
        }
        return dlist;
    }

    /**
     *
     * @param chr
     * @param pos
     * @return one fake {@link Digest} object (See {@link #makeFakeDigestList(String, Integer...)}).
     */
    private static Digest makeFakeDigest(String chr, int pos) {
        String fields[]=new String[6];
        fields[0]=chr;
        fields[1]=String.valueOf(pos);
        fields[2]=String.valueOf(pos+1);
        fields[3]=String.valueOf(++digestcounter);
        fields[4]=restrictionsite;
        fields[5]=restrictionsite;
       return new Digest(fields);
    }

    @Test
    public void testIteratorReturnsFirstPairOfReads() {
        sampairer = new SAMPairer(sam1,sam2,digestmap);
        Pair<SAMRecord,SAMRecord> pair = sampairer.getNextPair();
        assertNotNull(pair);
    }

    /**
     * XS:i:<n> Alignment score for the best-scoring alignment found other than the alignment reported.
     * Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode).
     * Only present if the SAM record is for an aligned read and more than one alignment was found for the read.
     * Note that, when the read is part of a concordantly-aligned pair, this score could be greater than AS:i.
     * The second read in reversetest.sam has these annotations: AS:i:0	XS:i:0, meaning that we want to filter out
     * this read. The first read should not be filtered out.
     */
    @Test
    public void testFilterMultipleAlignment() {
        sampairer = new SAMPairer(sam1,sam2,digestmap);
        Pair<SAMRecord,SAMRecord> pair = sampairer.getNextPair();
        assertTrue(sampairer.readPairUniquelyMapped(pair));
        pair=sampairer.getNextPair();
        assertFalse(sampairer.readPairUniquelyMapped(pair));
    }


    /**
     * An exception should be thrown if we look at a position that has no Digest.
     * @throws DiachromaticException
     */
    @Test(expected = DigestNotFoundException.class)
    public void testDigestNotFound() throws DiachromaticException{
        sampairer = new SAMPairer(sam1,sam2,digestmap);
        Pair<Digest, Digest> pair = sampairer.getDigestPair("crazyChromosome1", 1, 2, "wrongChromosome2", 3, 4);
    }
}
