package org.jax.diachromatic.map;

import com.sun.org.apache.regexp.internal.RE;
import htsjdk.samtools.SAMRecord;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.exception.DigestNotFoundException;
import org.jax.diachromatic.util.Pair;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.jax.diachromatic.map.TestSamFileImporter.retrieveReadPairsForTesting;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveSAMPairerForTesting;
import static org.junit.Assert.*;

public class SAMPairerTest {

    private Map<String,ReadPair> readpairmap;
    private SAMPairer sampairer;

    @Before
    public void init() {
        this.readpairmap = retrieveReadPairsForTesting();
        this.sampairer = retrieveSAMPairerForTesting();
    }

    /**
     * An exception should be thrown if we look at a position that has no Digest.
     * @throws DiachromaticException
     */
    @Test(expected = DigestNotFoundException.class)
    public void testDigestNotFound() throws DiachromaticException{
        DigestPair pair = sampairer.getDigestPair("crazyChromosome1", 1, 2, "wrongChromosome2", 3, 4);
    }


    @Test
    public void testDigestFound1() throws DiachromaticException {
        // we should find a digest for chr1	221618744 and for chr18	71915472
        // note the reads are 40bp long
        int READLEN=40;
        DigestPair pair = sampairer.getDigestPair("chr1",221618744,(221618744+READLEN),
                "chr18",71915472,(71915472 + READLEN));
        assertNotNull(pair);
    }

    @Test
    public void testDigestFound2() throws DiachromaticException {
        // we should find a digest for chr11 92316468 and for chr17	22262669
        int READLEN=40;
        //sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads);
        DigestPair pair = sampairer.getDigestPair("chr11",	92316468,(92316468+READLEN),
                "chr17",	22262669,(22262669+READLEN));
        assertNotNull(pair);
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
        ReadPair readpair =readpairmap.get("1_uniquelyAlignedRead");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        int insertSize=  readpair.getCalculatedInsertSize(digestpair);
        assertFalse(insertSize<LOWER_SIZE_THRESHOLD);
        assertFalse(insertSize>UPPER_SIZE_THRESHOLD);
        readpair = readpairmap.get("2_multiplyAlignedRead");
        digestpair = sampairer.getDigestPair(readpair);
        insertSize=  readpair.getCalculatedInsertSize(digestpair);
        assertFalse(insertSize<LOWER_SIZE_THRESHOLD);
        assertFalse(insertSize>UPPER_SIZE_THRESHOLD);
        readpair = readpairmap.get("3_tooBig");
        digestpair = sampairer.getDigestPair(readpair);
        insertSize=  readpair.getCalculatedInsertSize(digestpair);
        assertTrue(insertSize>UPPER_SIZE_THRESHOLD);
    }




    /** The sixth read pair is contiguous (by manual inspection) */
    @Test
    public void testContiguous() throws DiachromaticException {
        ReadPair readpair  = readpairmap.get("1_uniquelyAlignedRead");
        assertFalse(readpair.contiguous());
        readpair = readpairmap.get("5_religation"); //5 -- note readpair 5 was on adjacent fragments and
        // thus is religation and not contiguous!
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        assertTrue(readpair.religation(digestpair));
        //assertFalse(sampairer.contiguous(readpair))
        readpair = readpairmap.get("6_contiguous");//sampairer.getNextPair(); //6-- contiguous but not religated (not on adjacent digests)!
        assertTrue(readpair.contiguous());
        digestpair = sampairer.getDigestPair(readpair);
        assertFalse(readpair.religation(digestpair));
    }

    /** The insert of the third read pair is above threshold of 800. */
    @Test
    public void testInsertTooLarge() throws DiachromaticException {
        int THRESHOLD=800;
        ReadPair readpair  = readpairmap.get("3_tooBig");//1
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        int insertSize=readpair.getCalculatedInsertSize(digestpair);
        //System.err.println("insert size = " + insertSize); 3823
        assertTrue(insertSize>THRESHOLD);
    }

    /** The insert of the seventh read pair is above threshold of 800. */
    @Test
    public void testSetSamFlagsForCorrectPair() throws DiachromaticException {
        ReadPair readpair =readpairmap.get("7_validRead1");
        SamBitflagFilter.debugDisplayBitflag(readpair.forward().getFlags());
        // before we pair, the flags are set only to zero.
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
