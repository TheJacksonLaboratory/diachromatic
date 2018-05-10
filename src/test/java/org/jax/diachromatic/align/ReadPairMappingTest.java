package org.jax.diachromatic.align;


import org.jax.diachromatic.exception.DiachromaticException;
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

/**
 * This class is designed to test specifically the ability of the classes {@link Aligner} and
 * {@link ReadPair} to determine whether reads are uniquely mapped, unmapped, or multimapped.
 */
public class ReadPairMappingTest {
    private static int digestcounter=0;
    private static String restrictionsite=null;
    private static Map<String,List<Digest>> digestmap;
    private static Aligner sampairer;
    private static String sam1;
    private static String sam2;

    /** We will; use this to keep track of the read pairs from the test files forwardtest and reverse test.
     * Note that we have labeled the read pairs in those file to make them easy find for these tests.
     */
    private static Map<String,ReadPair> readpairmap;


    private static final boolean outputRejectedReads=false;

    @BeforeClass
    public  static void init() throws DiachromaticException {
        ClassLoader classLoader = MapPairsTest.class.getClassLoader();
        sam1 = classLoader.getResource("data/sam/forwardtest.sam").getFile();
        sam2 = classLoader.getResource("data/sam/reversetest.sam").getFile();
        // make fake digests for the reads we will test
        restrictionsite="HindIII";
        digestmap=new HashMap<>();/*
          String fields1[]= {"chr1","221612607","221618800","56554","HindIII","HindIII"};
        Digest d1 = new Digest(fields1);
        String fields2[]={"chr18","71910154","71919237","21136","HindIII",	"HindIII"};
        */
        List<Digest> digests = makeFakeDigestList("chr1",new Pair<>(221612607,221618800));
        digestmap.put("chr1",digests);
        // Note that we put in three contiguous segments, here. This is important for
        // test testContiguous(); if we leave out the middle segment, then the sixth read would
        // mistakenly be classified as religated (the reads are on the two outer Digests of a triplet!).
        digests = makeFakeDigestList("chr10", new Pair<>(53_177_087,53_179_904),
                new Pair<>(53_179_905,	53_180_587), new Pair<>(53_180_588,53187288));
        digestmap.put("chr10",digests);

        digests=makeFakeDigestList("chr11",new Pair<>(92314037,92316529),new Pair<>(92316530,92317574));
        digestmap.put("chr11",digests);
        // Note -- the following command will give these two digests neighboring digest numbers
        // this is needed for the test and resembles what we would have for a full dataset.
        digests=makeFakeDigestList("chr13",new Pair<>(31_421_583,31_425_191),new Pair<>(31425192,	31425873	));
        digestmap.put("chr13",digests);
        //chr15	91492555	91497580
        digests=makeFakeDigestList("chr15",new Pair<>(91_492_555,91_497_580));
        digestmap.put("chr15",digests);//chr15	91492809
        digests=makeFakeDigestList("chr16",new Pair<>(31497401,31527040), new Pair<>(84172259,84175274));
        digestmap.put("chr16",digests);
        digests=makeFakeDigestList("chr17",new Pair<>(22262265,22262874));
        digestmap.put("chr17",digests);
        digests=makeFakeDigestList("chr18",new Pair<>(71910154,71919237));
        digestmap.put("chr18",digests);
        readpairmap = new HashMap<>();
        String digestFile = "/Users/hansep/IdeaProjects/diachromatic/src/test/resources/data/testInteractionCountsMap/testInteractionCountsMapDigests.txt";
        String activeDigestsFile = "/Users/hansep/IdeaProjects/diachromatic/src/test/resources/data/testInteractionCountsMap/testInteractionCountsMapActiveDigests.txt";
        DigestMap digestMap = new DigestMap(digestFile, activeDigestsFile);
        sampairer = new Aligner(sam1,sam2,digestmap,outputRejectedReads,"test5",digestMap);
        ReadPair pair;
        while ((pair = sampairer.getNextPair())!=null) {
            readpairmap.put(pair.forward().getReadName(),pair);
            System.err.println(pair.forward().getReadName());
        }
    }

    /**
     * We require a digest list for the {@link Aligner} constructor. We make a "fake" digest list for
     * simplicity that will allow us to perform testing of various functionalities.
     * @param chrom
     * @param pos_pair
     * @return
     */
    private static List<Digest> makeFakeDigestList(String chrom,Pair<Integer,Integer> ...pos_pair) {
        List<Digest> dlist = new ArrayList<>();
        for (Pair<Integer,Integer> p : pos_pair) {
            Digest d = makeFakeDigest(chrom,p.first,p.second);
            dlist.add(d);
        }
        return dlist;
    }

    /**
     *
     * @param chr
     * @param topos
     * @param frompos
     * @return one fake {@link Digest} object (See {@link #makeFakeDigestList(String, Pair[])}).
     */
    private static Digest makeFakeDigest(String chr, int frompos, int topos) {
        String fields[]=new String[6];
        fields[0]=chr;
        fields[1]=String.valueOf(frompos);
        fields[2]=String.valueOf(topos);
        fields[3]=String.valueOf(++digestcounter);
        fields[4]=restrictionsite;
        fields[5]=restrictionsite;
        return new Digest(fields);
    }

    /**
     * We test the first read, {@code 1_uniquelyAlignedRead	0	chr16	31526917	8} and
     * {@code 1_uniquelyAlignedRead	0	chr16	84175204	42}. Neither of these reads have
     * the {@code XS:i:0} tag which would indicate multimapping.
     */
    @Test
    public void testUniquelyAlignedReadPair() {
        ReadPair pair=readpairmap.get("1_uniquelyAlignedRead");
        assertNotNull(pair);
        assertTrue(pair.readPairUniquelyMapped());
    }

    /**
     * The reads of the readpair {@code 2_multiplyAlignedRead} both have SAM bitfields of 0, but the
     * reverse read has a flag that indicates it is multimapped: {@code XS:i:0}.
     * <p>
     * XS:i:<n> Alignment score for the best-scoring alignment found other than the alignment reported.
     * Can be negative. Can be greater than 0 in --local mode (but not in --end-to-end mode).
     * Only present if the SAM record is for an aligned read and more than one alignment was found for the read.
     * Note that, when the read is part of a concordantly-aligned pair, this score could be greater than AS:i.
     * The second read in reversetest.sam has these annotations: AS:i:0	XS:i:0, meaning that we want to filter out
     * this read. The first read should not be filtered out.
     * </p>
     */
    @Test
    public void testMultiplyAlignedReadPair() {
        ReadPair pair=readpairmap.get("2_multiplyAlignedRead");
        assertNotNull(pair);
        assertFalse(pair.readPairUniquelyMapped());
        assertFalse(pair.isUnmapped());
        assertTrue(pair.isMultimapped());
    }


}
