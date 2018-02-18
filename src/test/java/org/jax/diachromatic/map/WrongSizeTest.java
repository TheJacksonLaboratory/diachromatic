package org.jax.diachromatic.map;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.Before;
import org.junit.Test;

import java.util.Map;
import java.util.Set;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveReadPairsForTesting;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveSAMPairerForTesting;

public class WrongSizeTest {

    private Map<String,ReadPair> readpairmap;
    private SAMPairer sampairer;
    int THRESHOLD=800;

    @Before
    public void init() {
        this.readpairmap = retrieveReadPairsForTesting();
        this.sampairer = retrieveSAMPairerForTesting();
    }

    /**
     * SRR071233.317.wrong_size_cis
     * forward: chr1:221618744
     * reverse: chr18:71915472
     * The insert size of this readpair is calculated to be 3823
     * @throws DiachromaticException
     */
    @Test
    public void testWrongSizeTrans() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.317.wrong_size_trans");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
        int insertSize=readpair.getCalculatedInsertSize(digestpair);
        //System.err.println("insert size = " + insertSize); //3823
        //
       // assertTrue(insertSize>THRESHOLD);
        assertFalse(readpair.hasValidInsertSize(digestpair));
        // also check we can adjust the threshold
        int originalThreshold = ReadPair.getUpperSizeThreshold();
        ReadPair.setUpperSizeThreshold(4000);
        assertTrue(readpair.hasValidInsertSize(digestpair));
        ReadPair.setUpperSizeThreshold(originalThreshold); // reset for next test!
    }

    /**
     * SRR071233.6637.wrong_size_cis
     * forward chr7	139942691 digest: chr7	139940523	139942837
     * reverse chr7	120599431 digest chr7	120595624	120605243
     * @throws DiachromaticException
     */
    @Test
    public void testWrongSizeCis() throws DiachromaticException {
        ReadPair readpair = readpairmap.get("SRR071233.6637.wrong_size_cis");
        assertNotNull(readpair);
        DigestPair digestpair = sampairer.getDigestPair(readpair);
//        int insertSize=readpair.getCalculatedInsertSize(digestpair);
//        System.err.println("insert size = " + insertSize); //3994
        assertFalse(readpair.hasValidInsertSize(digestpair));
    }






}
