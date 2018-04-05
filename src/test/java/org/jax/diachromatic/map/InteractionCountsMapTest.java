package org.jax.diachromatic.map;

import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;

import static org.junit.Assert.*;

public class InteractionCountsMapTest {

    private static InteractionCountsMap testInteractionCountsMap2c = null;
    private static InteractionCountsMap testInteractionCountsMap5c = null;

    @BeforeClass
    public static void setup() throws IOException {

        // create test objects for two conditions
        testInteractionCountsMap2c = new InteractionCountsMap(2);
        testInteractionCountsMap5c = new InteractionCountsMap(5);
    }

    @Test
    public void testUniqueKeys_IncrementFragPair() throws Exception {
        String uniqueKey=testInteractionCountsMap2c.incrementFragPair(0,"chr1",23, "chr3", 77);
        assertEquals(uniqueKey,"chr1:23:chr3:77");
        uniqueKey=testInteractionCountsMap2c.incrementFragPair(0,"chr3",23, "chr1", 77);
        assertEquals(uniqueKey,"chr3:23:chr1:77");
        uniqueKey=testInteractionCountsMap2c.incrementFragPair(0,"chr3",77, "chr1", 23);
        assertEquals(uniqueKey,"chr1:23:chr3:77");
        uniqueKey=testInteractionCountsMap2c.incrementFragPair(0,"chr1",77, "chr3", 23);
        assertEquals(uniqueKey,"chr3:23:chr1:77");
    }

    @Test
    public void testIncrementFragPair() throws Exception {
        String uniqueKey=testInteractionCountsMap5c.incrementFragPair(0,"chr1",23, "chr3", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(0,"chr1",23, "chr3", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(0,"chr1",23, "chr3", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(1,"chr3",23, "chr1", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(1,"chr3",23, "chr1", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(1,"chr3",23, "chr1", 77);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(2,"chr3",77, "chr1", 23);
        uniqueKey=testInteractionCountsMap5c.incrementFragPair(2,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.printInteractionCountsMap();
    }

    @Test
    public void testSameFragmentInteractionException() throws Exception {
        testInteractionCountsMap5c.incrementFragPair(2,"chr1",77, "chr1", 77);
        // TODO: Add proper exception handling and test.
    }

    @Test
    public void testGetCurrentTotalNumberOfInteractions() throws Exception {
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        assertEquals(6,testInteractionCountsMap5c.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        assertEquals(7,testInteractionCountsMap5c.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        assertEquals(8,testInteractionCountsMap5c.getCurrentTotalNumberOfInteractions().intValue());
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        assertEquals(9,testInteractionCountsMap5c.getCurrentTotalNumberOfInteractions().intValue());
    }

    @Test
    public void testGetInteractionNumForKeyAndCondition() {
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(0,"chr1",77, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(1,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.incrementFragPair(2,"chr2",99, "chr1", 23);
        testInteractionCountsMap5c.printInteractionCountsMap();
        System.out.println(testInteractionCountsMap5c.getInteractionNumForKeyAndCondition("chr1:23:chr1:77", 0));

    }
}