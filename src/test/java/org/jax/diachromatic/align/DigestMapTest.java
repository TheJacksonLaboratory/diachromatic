package org.jax.diachromatic.align;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import java.io.FileNotFoundException;
import java.net.URL;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

public class DigestMapTest {


    private static String digestFilePath;
    private static final double EPSILON=0.000001;

    @BeforeAll
    static void init() throws FileNotFoundException {
        ClassLoader classLoader = DigestMapTest.class.getClassLoader();
        URL url = classLoader.getResource("data/digestmap/hg19_digestedGenome.small.txt");
        if (url==null) {
            throw new FileNotFoundException("Could not find data/digestmap/hg19_digestGenome.small.txt");
        }
        digestFilePath = url.getFile();
    }

    @Test
    void testOpenFile() throws DiachromaticException  {
        DigestMap dmap = new DigestMap(digestFilePath);
        assertNotNull(dmap);
    }

    /**
     * hg19_digestedGenome.small.txt contains only chr1.
     * However, {@link DigestMap} adds both chr1 and 1
     */
    @Test
    void testGetOneChromosome() throws DiachromaticException {
        DigestMap dmap = new DigestMap(digestFilePath);
        Map<String, DigestMap.Chromosome2DigestArray> mymap = dmap.getDigestMap();
        int expectedNumberOfChromosomes = 2;
        assertEquals(expectedNumberOfChromosomes,mymap.size());
        assertTrue(mymap.containsKey("chr1"));
        assertTrue(mymap.containsKey("1"));
    }

    @Test
    void testGetDigests() throws DiachromaticException {
        DigestMap dmap = new DigestMap(digestFilePath);
        Map<String, DigestMap.Chromosome2DigestArray> mymap = dmap.getDigestMap();
        DigestMap.Chromosome2DigestArray chrom2array = mymap.get("chr1");
        int expectedNumberOfDigests = 19;
        Digest digest1 = chrom2array.getDigestAt(2);
        /*
        chr1	1	11159	1	None	DpnII	11159	0.000	0.000	0.000	0.022	F	0	0
         */
        assertEquals("chr1",digest1.getChromosome());
        assertEquals(1,digest1.getDigestStartPosition());
        assertEquals(11159,digest1.getDigestEndPosition());
        assertEquals("None",digest1.getFivePrimeRestrictionSite());
        assertEquals("DpnII",digest1.getThreePrimeRestrictionSite());
        assertEquals(0.000,digest1.getFivePrimeGcContent(),EPSILON);
        assertEquals(0.000,digest1.getThreePrimeGcContent(),EPSILON);
        assertEquals(0.000,digest1.getFivePrimeRepeatContent(),EPSILON);
        assertEquals(0.022,digest1.getThreePrimeRepeatContent(),EPSILON);
        assertFalse(digest1.isSelected());
        assertEquals(1,digest1.getDigesttNumber());
        assertEquals(0,digest1.getFivePrimeRepeatContent());
        assertEquals(0,digest1.getThreePrimeProbeCount());
    }


    // Test that the binary search gets the same result for any position in digest
    @Test
    void testGetSameDigest()throws DiachromaticException {
        //17954	18291
        DigestMap dmap = new DigestMap(digestFilePath);
        Map<String, DigestMap.Chromosome2DigestArray> mymap = dmap.getDigestMap();
        DigestMap.Chromosome2DigestArray chrom2array = mymap.get("chr1");

        Digest digest1 = chrom2array.getDigestAt(17954);
        Digest digest2 = chrom2array.getDigestAt(18000);
        Digest digest3 = chrom2array.getDigestAt(18291);
        assertEquals(digest1,digest2);
        assertEquals(digest1,digest3);
        assertEquals(digest2,digest3);
        // Not the same!
        Digest digest4 = chrom2array.getDigestAt(850);
        assertNotEquals(digest1,digest4);
    }

    @Test
    void testGetDigestPair()throws DiachromaticException {
        DigestMap dmap = new DigestMap(digestFilePath);
        Map<String, DigestMap.Chromosome2DigestArray> mymap = dmap.getDigestMap();
        DigestMap.Chromosome2DigestArray chrom2array = mymap.get("chr1");
        Digest digest1 = chrom2array.getDigestAt(17954);
        Digest digest4 = chrom2array.getDigestAt(850);

        DigestPair dpair = dmap.getDigestPair("chr1",17954,"chr1",850);
        assertEquals(digest1,dpair.forward());
        assertEquals(digest4,dpair.reverse());
    }

}
