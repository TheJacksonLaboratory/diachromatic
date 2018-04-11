package org.jax.diachromatic.map;

import org.jax.diachromatic.exception.IncrementSameInternalInteractionException;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.io.IOException;

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
    public void testUniqueKeys_IncrementFragPair() throws IncrementSameInternalInteractionException {
        String uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 23, 33,"chr3", 77,87);
        assertEquals(uniqueKey, "chr1:23-33;chr3:77-87");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 23, 33,"chr1", 77, 87);
        assertEquals(uniqueKey, "chr3:23-33;chr1:77-87");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 77, 87,"chr1", 23, 33);
        assertEquals(uniqueKey, "chr1:23-33;chr3:77-87");
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 77, 87, "chr3", 23, 33);
        assertEquals(uniqueKey, "chr3:23-33;chr1:77-87");
    }

    /**
     * Test if fields are incremented as expected.
     *
     * @throws Exception
     */
    @Test
    public void testIncrementFragPair() throws IncrementSameInternalInteractionException, FileNotFoundException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);
        String uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23, 33,"chr3", 77, 87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23,33,"chr3", 77,87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23,33,"chr3", 77,87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(2, "chr3", 77,87,"chr1", 23,33);
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(2, "chr1", 77,87,"chr1", 23,33);
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr1:77-87", 0).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr1:77-87", 1).intValue());
        assertEquals(1, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr1:77-87", 2).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr1:77-87", 3).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr1:77-87", 4).intValue());

        assertEquals(3, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr3:77-87", 0).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr3:77-87", 1).intValue());
        assertEquals(1, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr3:77-87", 2).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr3:77-87", 3).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr1:23-33;chr3:77-87", 4).intValue());

        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr3:23-33;chr1:77-87", 0).intValue());
        assertEquals(3, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr3:23-33;chr1:77-87", 1).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr3:23-33;chr1:77-87", 2).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr3:23-33;chr1:77-87", 3).intValue());
        assertEquals(0, testInteractionCountsMap5c.getNumberOfInteractionsForKeyAndCondition("chr3:23-33;chr1:77-87", 4).intValue());
        testInteractionCountsMap5c.printInteractionCountsMapAsCountTable();
    }

    /**
     * Test exception for the incrementation of interaction counts between identical fragments.
     *
     * @throws Exception
     *
     */
    @Test
    public void testSameFragmentInteractionException() throws IncrementSameInternalInteractionException {
        testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77, 87,"chr1", 77, 87);
    }

    /**
     * Test if total number of interactions is counted correctly.
     *
     * @throws Exception
     */
    @Test
    public void testGetTotalNumberOfInteractionsForCondition() throws IncrementSameInternalInteractionException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,"chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,"chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,"chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109,"chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109, "chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109, "chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, 109,"chr1", 23, 33);
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, 109,"chr1", 23, 33);
        assertEquals(3, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(0).intValue());
        assertEquals(3, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(1).intValue());
        assertEquals(2, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(2).intValue());
    }

    @Test
    public void printInteractionCountsMapCountTable() throws IncrementSameInternalInteractionException, FileNotFoundException {
        String uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23, 33,"chr3", 77, 87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23,33,"chr3", 77,87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23,33,"chr3", 77,87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,"chr1", 77,87);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr3", 77,87,"chr1", 23,33);
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77,87,"chr1", 23,33);
        testInteractionCountsMap5c.printInteractionCountsMapAsCountTable();
    }

    /**
     * Test if the numbers of reads at interacting fragments are derived correctly from interaction counts.
     *
     * @throws FileNotFoundException
     * @throws IncrementSameInternalInteractionException
     */
    @Test
    public void testDeriveReadCountsAtInteractingFragments() throws FileNotFoundException, IncrementSameInternalInteractionException {

        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);

        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,"chr1", 40, 50);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,"chr1", 70, 80);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,"chr1", 90, 100);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,"chr1", 110, 120);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,"chr1", 110, 120);
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,"chr1", 110, 120);

        testInteractionCountsMapTmp.deriveReadCountsAtInteractingFragments();

        assertEquals(4,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:90-100",0).intValue());
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:40-50",0).intValue());
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:70-80",0).intValue());
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:110-120",0).intValue());
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:10-20",0).intValue());

        testInteractionCountsMapTmp.printFragmentInteractionCountsMapAsCountTable();
        testInteractionCountsMapTmp.printInteractionCountsMapAsCountTable();

        System.out.println("TotalNumberOfInteractionsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(0));
        System.out.println("TotalNumberOfInteractingFragmentsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractingFragmentsForCondition(0));
    }
}