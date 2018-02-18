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
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

public class ReligationAndContiguousTest {
    private Map<String,ReadPair> readpairmap;
    private SAMPairer sampairer;

    @Before
    public void init() {
        this.readpairmap = retrieveReadPairsForTesting();
        this.sampairer = retrieveSAMPairerForTesting();
    }

    /** The read pair SRR071233.18645.re_ligation is not on the same fragment.
     * Instead, the reads are on adjacent fragments and have religated.
     * @throws DiachromaticException
     */
    @Test
    public void testReligation() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.18645.re_ligation");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        boolean sameRestrictionFragment=readpair.bothReadsLocatedOnSameRestrictionFragment(digestpair);
        assertFalse(sameRestrictionFragment);
        assertTrue(readpair.religation(digestpair));
    }


    /* The fifth pair shows religation! */
    @Test
    public void testReligation2() throws DiachromaticException {

        ReadPair readpair  = readpairmap.get("1_uniquelyAlignedRead");
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        assertFalse(readpair.religation(digestpair));
        readpair = readpairmap.get("4_selfLigation");// fourth read pair, self-ligation--not religation
        digestpair = sampairer.getDigestPair(readpair);
        assertFalse(readpair.religation(digestpair));
        //chr13	31421583	31425191
        //chr13	31425192	31425873
        readpair = readpairmap.get("5_religation");// 5. religation!
        assertEquals(readpair.forward().getReferenceName(),"chr13");// check we have right read!
        assertEquals(readpair.reverse().getReferenceName(),"chr13");// check we have right read!
        digestpair = sampairer.getDigestPair(readpair);
        Assert.assertTrue(readpair.religation(digestpair));
    }




}
