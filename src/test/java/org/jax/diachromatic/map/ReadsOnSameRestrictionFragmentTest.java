package org.jax.diachromatic.map;

import org.jax.diachromatic.exception.DiachromaticException;
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
        Set<QCCode> qcCodes=readpair.getErrorCodes();
        assertTrue(qcCodes.contains(QCCode.SAME_DANGLING_END));
    }


    /**
     * SRR071233.19593.same_circularised	has these reads
     * forward: chr15:91_492_809 TGCTGAGCCTAGGAGTGCAGCGTCTGCCTGCCAGAGTCCC
     * reverse: chr15:91_497_543 AGCTAATCAATAAGAAACTCACCAGGTGCGGTGCCTCAAG
     * Both reads are on the same digest fragment chr15:91_492_555-91_497_580
     * Note that the reverse fragment overlaps the end of the digest fragment (it extends two nucleotides
     * beyond 555-91_497_580, as can be seen in the UCSC browser).
     * @throws DiachromaticException
     */
    @Test
    public void testSameCircularized1() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.19593.same_circularised");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        boolean sameRestrictionFragment=readpair.bothReadsLocatedOnSameRestrictionFragment(digestpair);
        assertTrue(sameRestrictionFragment);
        Set<QCCode> qcCodes=readpair.getErrorCodes();

        assertTrue(qcCodes.contains(QCCode.CIRCULARIZATION_DANGLING));
    }

    /**
     * SRR071233.24017.dangling_ends has
     * forward: chr7	34_167_256 (SAM code 16=reverse strand)
     * reverse: chr7	34_167_161 (SAM code 0=forward strand)
     * Both are on the same digest: chr7:34_161_770-34_167_291
     * These reads both "point to the middle" and thus are not circularized. the forward read overlaps with the
     * digest. They are therefore same dangling end
     * @throws DiachromaticException
     */
    @Test
    public void testDangling() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.24017.dangling_ends");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        boolean sameRestrictionFragment=readpair.bothReadsLocatedOnSameRestrictionFragment(digestpair);
        assertTrue(sameRestrictionFragment);
        Set<QCCode> qcCodes=readpair.getErrorCodes();
        for (QCCode q : qcCodes) {
            System.err.println(q);
        }
        assertTrue(qcCodes.contains(QCCode.SAME_DANGLING_END));
        assertTrue(readpair.hasValidInsertSize(digestpair));
    }

    /**
     * SRR071233.2213.same_internal
     * forward: chr6:67_518_942  (SAM code 0=forward strand)
     * reverse: chr6: 67_519_165 (SAM code 16=reverse strand)
     * Both have the same digest, chr6:67_509_031-67_520_069
     * These reads both "point to the middle" and thus are not circularized. They do not overlap with the
     * digest ends. They are therefore same internal
     * @throws DiachromaticException
     */
    @Test
    public void testSameInternal() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.2213.same_internal");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        boolean sameRestrictionFragment=readpair.bothReadsLocatedOnSameRestrictionFragment(digestpair);
        assertTrue(sameRestrictionFragment);
        Set<QCCode> qcCodes=readpair.getErrorCodes();
        assertTrue(qcCodes.contains(QCCode.SAME_INTERNAL));
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
