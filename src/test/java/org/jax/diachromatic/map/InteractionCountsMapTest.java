package org.jax.diachromatic.map;

import org.jax.diachromatic.exception.IncrementSameInternalInteraction;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;
import static org.junit.Assert.assertEquals;

public class InteractionCountsMapTest {

    private static InteractionCountsMap testInteractionCountsMap2c = null;
    private static InteractionCountsMap testInteractionCountsMap5c = null;

    @BeforeClass
    public static void setup() throws IOException {

        // create test objects for two conditions
        testInteractionCountsMap2c = new InteractionCountsMap(2);
        testInteractionCountsMap5c = new InteractionCountsMap(5);
    }

    /**
     * Test if keys are properly assembled. The fragment with the smaller coordinate comes always first.
     *
     * @throws Exception
     */
    @Test
    public void testUniqueKeys_IncrementFragPair() throws Exception {
        String uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 23, "chr3", 77);
        assertEquals(uniqueKey, "chr1:23:chr3:77");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 23, "chr1", 77);
        assertEquals(uniqueKey, "chr3:23:chr1:77");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 77, "chr1", 23);
        assertEquals(uniqueKey, "chr1:23:chr3:77");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 77, "chr3", 23);
        assertEquals(uniqueKey, "chr3:23:chr1:77");
    }

    /**
     * Test if fields are incremented as expected.
     *
     * @throws Exception
     */
    @Test
    public void testIncrementFragPair() throws Exception {
        String uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23, "chr3", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23, "chr3", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23, "chr3", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23, "chr1", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23, "chr1", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23, "chr1", 77);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr3", 77, "chr1", 23);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77, "chr1", 23);
        //testInteractionCountsMap5c.printInteractionCountsMap();
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 0).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 1).intValue());
        assertEquals(1, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 2).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 3).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 4).intValue());

        assertEquals(3, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr3:77", 0).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr3:77", 1).intValue());
        assertEquals(1, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr3:77", 2).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr3:77", 3).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr3:77", 4).intValue());

        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr3:23:chr1:77", 0).intValue());
        assertEquals(3, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr3:23:chr1:77", 1).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr3:23:chr1:77", 2).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr3:23:chr1:77", 3).intValue());
        assertEquals(0, testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr3:23:chr1:77", 4).intValue());
    }

    /**
     * Test exception for the incrementation of interaction counts between identical fragments.
     *
     * @throws Exception
     *
     * TODO: Find a way in JUnit to test if IncrementSameInternalInteraction is thrown.
     *
     */
    @Test
    public void testSameFragmentInteractionException() throws Exception {
        testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77, "chr1", 77);
    }

    /**
     * Test if total number of interactions is counted correctly.
     *
     * @throws Exception
     */
    @Test
    public void testGetCurrentTotalNumberOfInteractions() throws Exception {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, "chr1", 23);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, "chr1", 23);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, "chr1", 23);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, "chr1", 23);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, "chr1", 23);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, "chr1", 23);
        assertEquals(6, testInteractionCountsMapTmp.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, "chr1", 23);
        assertEquals(7, testInteractionCountsMapTmp.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, "chr1", 23);
        assertEquals(8, testInteractionCountsMapTmp.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, "chr1", 23);
        assertEquals(9, testInteractionCountsMapTmp.getCurrentTotalNumberOfInteractions().intValue());
    }
}