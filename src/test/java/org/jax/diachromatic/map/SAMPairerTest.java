package org.jax.diachromatic.map;

import com.sun.org.apache.regexp.internal.RE;
import htsjdk.samtools.SAMRecord;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;
import org.jax.diachromatic.util.Pair;
import org.junit.BeforeClass;
import org.junit.Test;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.*;

public class SAMPairerTest {

    private static int digestcounter=0;
    private static String restrictionsite=null;
    private static Map<String,List<Digest>> digestmap;
    private static SAMPairer sampairer;
    private static String sam1;
    private static String sam2;

    /** We will; use this to keep track of the read pairs from the test files forward test and reverse test.
     * Note that we have labeled the read pairs in those file to make them easy find for these tests.
     */
    private static Map<String,ReadPair> readpairmap;


    private static final boolean outputRejectedReads=false;

    @BeforeClass
    public  static void init() throws DiachromaticException {
        ClassLoader classLoader = SAMPairerTest.class.getClassLoader();
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
        String outdir = "results";
        String outprefix = "results";
        sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads,outdir, outprefix);
        ReadPair pair;
        while ((pair = sampairer.getNextPair())!=null) {
            readpairmap.put(pair.forward().getReadName(),pair);
            System.err.println(pair.forward().getReadName());
        }
    }


    /**
     * We require a digest list for the {@link SAMPairer} constructor. We make a "fake" digest list for
     * simplicity that will allow us to perform testing of various functionalities.
     * @param chrom chromosome
     * @param pos_pair pairs of integers with start and end position of the Digests
     * @return list of "fake" Digest object
     */
    @SafeVarargs
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
     * @param chr chromosome
     * @param topos start position of digest
     * @param frompos end position of digest
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
     * We are testing the third pair of reads, which was extracted from the file
     * {@code _wrong_size.filter.sam}. The first two
     * reads should be OK. The digests are created to encompass these two reads
     * from the corresponding digest file.
     * We are testing the fragments are ok in size. Note that we have added the corresponding
     * Digests for the first three fragments.
     */
    @Test
    public void testFragmentToLarge() throws DiachromaticException {
        int UPPER_SIZE_THRESHOLD=800;
        int LOWER_SIZE_THRESHOLD=150;
        sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads,"results", "prefix");
        ReadPair readpair =readpairmap.get("1_uniquelyAlignedRead");
        assertNotNull(readpair);
        int insertSize=  readpair.getCalculatedInsertSize();
        assertFalse(insertSize<LOWER_SIZE_THRESHOLD);
        assertFalse(insertSize>UPPER_SIZE_THRESHOLD);
        readpair = readpairmap.get("3_tooBig");
        insertSize=  readpair.getCalculatedInsertSize();
        assertTrue(insertSize>UPPER_SIZE_THRESHOLD);
    }

    /** This should throw an exception because we cannot calculate the insert size for a multimapped read. */
    @Test(expected =DiachromaticException.class)
    public void testFragmentToLargeException() throws DiachromaticException {
        int UPPER_SIZE_THRESHOLD=800;
        int LOWER_SIZE_THRESHOLD=150;
        sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads,"results", "prefix");
        ReadPair readpair = readpairmap.get("2_multiplyAlignedRead");
        assertNotNull(readpair);
        int insertSize=  readpair.getCalculatedInsertSize();
        assertFalse(insertSize<LOWER_SIZE_THRESHOLD);
        assertFalse(insertSize>UPPER_SIZE_THRESHOLD);
    }

    /**
     * We are testing the fourth pair of reads that self-circularizes. The first read
     * pair is from the same chromosome but does not self-circularize. (by manual inspection).
     * The second and third pairs have distinct chromosomes and cnnot be tested for self circularization.
     * @throws DiachromaticException
     */
    @Test
    public void testSelfLigation() throws DiachromaticException {
        sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads,"results", "prefix");
        ReadPair readpair = readpairmap.get("1_uniquelyAlignedRead");
        assertFalse(readpair.selfLigation());
        readpair = readpairmap.get("4_selfLigation");//sampairer.getNextPair();// fourth read pair, self-ligation!
        assertTrue(readpair.selfLigation());
    }

    /* The fifth pair shows religation! */
    @Test
    public void testReligation() throws DiachromaticException {

        ReadPair readpair  = readpairmap.get("1_uniquelyAlignedRead");
        assertFalse(readpair.religation());
        readpair = readpairmap.get("4_selfLigation");// fourth read pair, self-ligation--not religation
        assertFalse(readpair.religation());
        //chr13	31421583	31425191
        //chr13	31425192	31425873
        readpair = readpairmap.get("5_religation");// 5. religation!
        assertEquals(readpair.forward().getReferenceName(),"chr13");// check we have right read!
        assertEquals(readpair.reverse().getReferenceName(),"chr13");// check we have right read!
        assertTrue(readpair.religation());
    }

    /** The sixth read pair is contiguous (by manual inspection) */
    @Test
    public void testIsContiguous() throws DiachromaticException {
        ReadPair readpair  = readpairmap.get("1_uniquelyAlignedRead");
        assertFalse(readpair.isContiguous());
        readpair = readpairmap.get("5_religation"); //5 -- note readpair 5 was on adjacent fragments and
        // thus is religation and not contiguous!
        assertTrue(readpair.religation());
        //assertFalse(sampairer.contiguous(readpair))
        readpair = readpairmap.get("6_contiguous");//sampairer.getNextPair(); //6-- contiguous but not religated (not on adjacent digests)!
        assertTrue(readpair.isContiguous());
        assertFalse(readpair.religation());
    }

    /** The insert of the third read pair is above threshold of 800. */
    @Test
    public void testInsertTooLarge() throws DiachromaticException {
        int THRESHOLD=800;
        ReadPair readpair  = readpairmap.get("3_tooBig");//1
        int insertSize=readpair.getCalculatedInsertSize();
        //System.err.println("insert size = " + insertSize); 3823
        assertTrue(insertSize>THRESHOLD);
    }


    /** The insert of the seventh read pair is above threshold of 800. */
    @Test
    public void testSetSamFlagsForCorrectPair() throws DiachromaticException {
        ReadPair readpair =readpairmap.get("7_validRead1");
        SamBitflagFilter.debugDisplayBitflag(readpair.forward().getFlags());
        // before we pair, the flags are set only to zero.
        readpair.forward().setFlags(0);
        readpair.reverse().setFlags(0);
        assertEquals(0,readpair.forward().getFlags());
        assertEquals(0,readpair.reverse().getFlags());
        readpair.pairReads();
        assertEquals(67,readpair.forward().getFlags());
        assertEquals(131,readpair.reverse().getFlags());
        SamBitflagFilter.debugDisplayBitflag(readpair.forward().getFlags());
    }


    //TODO CHECK THIS TEST
    @Test
    public void testDuplicate() {
        ReadPair readpair =readpairmap.get("7_validRead1");
        assertFalse(DiTag.isDuplicate(readpair));
        readpair=readpairmap.get("8_validRead1");
//        assertTrue(DiTag.isDuplicate(readpair));
    }



}
