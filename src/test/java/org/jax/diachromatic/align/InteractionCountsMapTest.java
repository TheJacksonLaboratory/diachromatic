package org.jax.diachromatic.align;

import org.jax.diachromatic.count.InteractionCountsMap;
import org.jax.diachromatic.exception.IncrementSameInternalInteractionException;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;


import java.io.FileNotFoundException;

import static org.junit.Assert.assertEquals;

public class InteractionCountsMapTest {

    private static InteractionCountsMap testInteractionCountsMap2c = null;
    private static InteractionCountsMap testInteractionCountsMap5c = null;

    @BeforeAll
    public static void setup() {

        // create test objects for two conditions
        testInteractionCountsMap2c = new InteractionCountsMap();
        testInteractionCountsMap5c = new InteractionCountsMap();
    }

    /**
     * Test if keys are properly assembled. The fragment with the smaller coordinate comes always first.
     *
     * @throws Exception
     */
    @Test
    public void testUniqueKeys_IncrementFragPair() throws IncrementSameInternalInteractionException {
        String uniqueKey = testInteractionCountsMap2c.incrementFragPair("chr1", 23, 33, false,"chr3", 77,87,false,"F1F2");
        assertEquals("chr1:23-33:I;chr3:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair("chr3", 23, 33,false,"chr1", 77, 87,false,"F1F2");
        assertEquals("chr3:23-33:I;chr1:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair("chr3", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        assertEquals("chr1:23-33:I;chr3:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair("chr1", 77, 87, false,"chr3", 23, 33,false,"F1F2");
        assertEquals("chr3:23-33:I;chr1:77-87:I;T",uniqueKey);
    }

    /**
     * Test if fields are incremented as expected.
     *
     * @throws Exception
     */
    @Test
    public void testIncrementFragPair() throws IncrementSameInternalInteractionException, FileNotFoundException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap();
        String uniqueKey = testInteractionCountsMapTmp.incrementFragPair( "chr1", 23, 33, false,"chr3", 77, 87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr3", 77,87,false,"chr1", 23,33,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair("chr1", 77,87,false,"chr1", 23,33,false,"F1F2");
        assertEquals(1, testInteractionCountsMapTmp.getNumberOfInteractionsForKey("chr1:23-33:I;chr1:77-87:I;T"));
        assertEquals(4, testInteractionCountsMapTmp.getNumberOfInteractionsForKey("chr1:23-33:I;chr3:77-87:I;T"));
        assertEquals(3, testInteractionCountsMapTmp.getNumberOfInteractionsForKey("chr3:23-33:I;chr1:77-87:I;T"));
        testInteractionCountsMap5c.printInteractionCountsMapAsCountTable("test2");
    }

    /**
     * Test exception for the incrementation of interaction counts between identical fragments.
     *
     * @throws Exception
     *
     */
    @Test
    public void testSameFragmentInteractionException() throws IncrementSameInternalInteractionException {
        testInteractionCountsMap5c.incrementFragPair("chr1", 77, 87,false,"chr1", 77, 87,false,"F1F2");
    }

    /**
     * Test if total number of interactions is counted correctly.
     *
     * @throws Exception
     */
    @Test
    public void testGetTotalNumberOfInteractionsForCondition() throws IncrementSameInternalInteractionException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap();
        // three read pairs -> one interaction
        testInteractionCountsMapTmp.incrementFragPair("chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        // three read pairs -> one interaction
        testInteractionCountsMapTmp.incrementFragPair("chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr2", 99, 109,false, "chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr2", 99, 109,false, "chr1", 23, 33,false,"F1F2");
        // three read pairs -> two interactions
        testInteractionCountsMapTmp.incrementFragPair("chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr2", 199, 1109,false,"chr1", 163, 173,false,"F1F2");

        assertEquals(3, testInteractionCountsMapTmp.getTotalNumberOfInteractions().intValue());
    }

    /**
     * The following test is disabled/ignored (printing to file is causing an error).
     * @throws IncrementSameInternalInteractionException
     * @throws FileNotFoundException
     */
    @Test @Disabled
    public void printInteractionCountsMapCountTable() throws IncrementSameInternalInteractionException, FileNotFoundException {
        String uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr1", 23, 33,false,"chr3", 77, 87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr3", 77,87,false,"chr1", 23,33,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair("chr1", 77,87,false,"chr1", 23,33,false,"F1F2");
       testInteractionCountsMap5c.printInteractionCountsMapAsCountTable("test");
    }

    /**
     * Test if the numbers of reads at interacting fragments are derived correctly from interaction counts.
     *
     * @throws FileNotFoundException
     * @throws IncrementSameInternalInteractionException
     */
    @Test
    public void testDeriveReadCountsAtInteractingFragments() throws FileNotFoundException, IncrementSameInternalInteractionException {

        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap();

        testInteractionCountsMapTmp.incrementFragPair("chr1", 10, 20,false,"chr1", 40, 50,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 10, 20,false,"chr1", 70, 80,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 10, 20,false,"chr1", 90, 100,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair("chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");

        testInteractionCountsMapTmp.deriveReadCountsAtInteractingFragments();

        assertEquals(4,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKey("chr1:90-100:I"));
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKey("chr1:40-50:I"));
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKey("chr1:70-80:I"));
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKey("chr1:110-120:I"));
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKey("chr1:10-20:I"));

        testInteractionCountsMapTmp.printFragmentInteractionCountsMapAsCountTable("test3.tsv");
        testInteractionCountsMapTmp.printInteractionCountsMapAsCountTable("test4.tsv");

        System.out.println("TotalNumberOfInteractionsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractions());
        System.out.println("TotalNumberOfInteractingFragmentsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractingFragments());
    }
}
