package org.jax.diachromatic.map;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Map;
import java.util.Set;

import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveReadPairsForTesting;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveSAMPairerForTesting;
import static org.junit.Assert.assertFalse;

/** This test class is intended to test whether we can recognize dangling ends.
 * It is also possible for non-ligated DNA fragments to insert between sequencing
 * adapters, despite the protocol being designed to minimise such events. Consequently,
 * HiCUP identifies and removes these unwanted species by checking if the forward and
 * reverse reads of a di-tag map to the same genomic restriction fragment, but unlike
 * circularised fragments the reads are orientated towards each other. Furthermore,
 * HiCUP divides this category into two sub-categories depending on whether the DNA fragment
 * end overlaps a restriction fragment cut site.
 * If the fragment end does overlap it is termed a "dangling end" (Figure 2d),
 * but if it does not it is termed an "internal fragment" (Figure 2e).
 */
public class ReadsOnSameRestrictionFragmentTest {

    private Map<String,ReadPair> readpairmap;
    private SAMPairer sampairer;

    @Before
    public void init() {
        this.readpairmap = retrieveReadPairsForTesting();
        this.sampairer = retrieveSAMPairerForTesting();
    }


    @Test
    public void testDanglingEnds() throws DiachromaticException{
        ReadPair readpair = readpairmap.get("SRR071233.24017.dangling_ends");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        boolean sameRestrictionFragment=readpair.bothReadsLocatedOnSameRestrictionFragment(digestpair);
        assertTrue(sameRestrictionFragment);
        Set<ErrorCode> qcCodes=readpair.getErrorCodes();
        assertTrue(qcCodes.contains(ErrorCode.SAME_DANGLING_END));
    }

    /**
     * We are testing the fourth pair of reads that self-circularizes. The first read
     * pair is from the same chromosome but does not self-circularize. (by manual inspection).
     * The second and third pairs have distinct chromosomes and cnnot be tested for self circularization.
     * @throws DiachromaticException
     */
//    @Test
//    public void testSelfLigation() throws DiachromaticException {
//        sampairer = new SAMPairer(sam1,sam2,digestmap,outputRejectedReads);
//        ReadPair readpair = readpairmap.get("1_uniquelyAlignedRead");
//        assertFalse(readpair.selfLigation());
//        readpair = readpairmap.get("4_selfLigation");//sampairer.getNextPair();// fourth read pair, self-ligation!
//        Assert.assertTrue(readpair.selfLigation());
//    }


}
