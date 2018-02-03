package org.jax.diachromatic.map;

import org.junit.Before;
import org.junit.Test;

import java.util.Map;

import static junit.framework.TestCase.assertNotNull;
import static org.jax.diachromatic.map.TestSamFileImporter.retrieveReadPairsForTesting;

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
public class DanglingEndTest {

    private Map<String,ReadPair> readpairmap;

    @Before
    public void init() {
        this.readpairmap = retrieveReadPairsForTesting();
    }


    @Test
    public void testDanglingEnds() {
        ReadPair readpair = readpairmap.get("SRR071233.24017.dangling_ends");
        assertNotNull(readpair);
    }


}
