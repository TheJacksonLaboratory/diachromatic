package org.jax.diachromatic.align;

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
        String uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 23, 33, false,"chr3", 77,87,false,"F1F2");
        assertEquals("chr1:23-33:I;chr3:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 23, 33,false,"chr1", 77, 87,false,"F1F2");
        assertEquals("chr3:23-33:I;chr1:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr3", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        assertEquals("chr1:23-33:I;chr3:77-87:I;T",uniqueKey);
        uniqueKey = testInteractionCountsMap2c.incrementFragPair(0, "chr1", 77, 87, false,"chr3", 23, 33,false,"F1F2");
        assertEquals("chr3:23-33:I;chr1:77-87:I;T",uniqueKey);
    }

    /**
     * Test if fields are incremented as expected.
     *
     * @throws Exception
     */
    @Test
    public void testIncrementFragPair() throws IncrementSameInternalInteractionException, FileNotFoundException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);
        String uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23, 33, false,"chr3", 77, 87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(2, "chr3", 77,87,false,"chr1", 23,33,false,"F1F2");
        uniqueKey = testInteractionCountsMapTmp.incrementFragPair(2, "chr1", 77,87,false,"chr1", 23,33,false,"F1F2");
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr1:77-87:I;T", 0).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr1:77-87:I;T", 1).intValue());
        assertEquals(1, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr1:77-87:I;T", 2).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr1:77-87:I;T", 3).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr1:77-87:I;T", 4).intValue());

        assertEquals(3, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr3:77-87:I;T", 0).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr3:77-87:I;T", 1).intValue());
        assertEquals(1, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr3:77-87:I;T", 2).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr3:77-87:I;T", 3).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr1:23-33:I;chr3:77-87:I;T", 4).intValue());

        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr3:23-33:I;chr1:77-87:I;T", 0).intValue());
        assertEquals(3, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr3:23-33:I;chr1:77-87:I;T", 1).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr3:23-33:I;chr1:77-87:I;T", 2).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr3:23-33:I;chr1:77-87:I;T", 3).intValue());
        assertEquals(0, testInteractionCountsMapTmp.getNumberOfInteractionsForKeyAndCondition("chr3:23-33:I;chr1:77-87:I;T", 4).intValue());
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
        testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77, 87,false,"chr1", 77, 87,false,"F1F2");
    }

    /**
     * Test if total number of interactions is counted correctly.
     *
     * @throws Exception
     */
    @Test
    public void testGetTotalNumberOfInteractionsForCondition() throws IncrementSameInternalInteractionException {
        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);
        // three read pairs -> one interaction
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 77, 87,false,"chr1", 23, 33,false,"F1F2");
        // three read pairs -> one interaction
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109,false, "chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(1, "chr2", 99, 109,false, "chr1", 23, 33,false,"F1F2");
        // three read pairs -> two interactions
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 99, 109,false,"chr1", 23, 33,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(2, "chr2", 199, 1109,false,"chr1", 163, 173,false,"F1F2");

        assertEquals(1, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(0).intValue());
        assertEquals(1, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(1).intValue());
        assertEquals(2, testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(2).intValue());
    }

    @Test
    public void printInteractionCountsMapCountTable() throws IncrementSameInternalInteractionException, FileNotFoundException {
        String uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23, 33,false,"chr3", 77, 87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(0, "chr1", 23,33,false,"chr3", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(1, "chr3", 23,33,false,"chr1", 77,87,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr3", 77,87,false,"chr1", 23,33,false,"F1F2");
        uniqueKey = testInteractionCountsMap5c.incrementFragPair(2, "chr1", 77,87,false,"chr1", 23,33,false,"F1F2");
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

        InteractionCountsMap testInteractionCountsMapTmp = new InteractionCountsMap(5);

        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,false,"chr1", 40, 50,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,false,"chr1", 70, 80,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 10, 20,false,"chr1", 90, 100,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");
        testInteractionCountsMapTmp.incrementFragPair(0, "chr1", 90, 100,false,"chr1", 110, 120,false,"F1F2");

        testInteractionCountsMapTmp.deriveReadCountsAtInteractingFragments();

        assertEquals(4,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:90-100:I",0).intValue());
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:40-50:I",0).intValue());
        assertEquals(1,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:70-80:I",0).intValue());
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:110-120:I",0).intValue());
        assertEquals(3,testInteractionCountsMapTmp.getNumberOfReadsAtInteractingFragmentForKeyAndCondition("chr1:10-20:I",0).intValue());

        testInteractionCountsMapTmp.printFragmentInteractionCountsMapAsCountTable("test3.tsv");
        testInteractionCountsMapTmp.printInteractionCountsMapAsCountTable("test4.tsv");

        System.out.println("TotalNumberOfInteractionsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractionsForCondition(0));
        System.out.println("TotalNumberOfInteractingFragmentsForCondition: " + testInteractionCountsMapTmp.getTotalNumberOfInteractingFragmentsForCondition(0));
    }
}