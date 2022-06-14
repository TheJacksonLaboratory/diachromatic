package org.jax.diachromatic.align;


import htsjdk.samtools.SAMRecord;
import org.junit.Ignore;
import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;


import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;



@Disabled("Cannot run Mockito with final class from HTSJDK")
class DeDupMapTest {


    @Test
    void recognizeDuplicateTest() {
        boolean useRelativeOrientation = false;
        DeDupMap ddmap = new DeDupMap(useRelativeOrientation);
        ReadPair rpair1 = Mockito.mock(ReadPair.class);
        // mock SAMRecords, forward and reverse for the two read pairs.
        SAMRecord r1f = Mockito.mock(SAMRecord.class);
        SAMRecord r1r = Mockito.mock(SAMRecord.class);
        SAMRecord r2f = Mockito.mock(SAMRecord.class);
        SAMRecord r2r = Mockito.mock(SAMRecord.class);
        when(r1f.getReferenceName()).thenReturn("chrZ");
        when(r1r.getReferenceName()).thenReturn("chrZ");
        when(r2f.getReferenceName()).thenReturn("chrZ");
        when(r2r.getReferenceName()).thenReturn("chrZ");

        when(rpair1.isTrans()).thenReturn(false);
        when(rpair1.getFivePrimeEndPosOfR1()).thenReturn(10);
        when(rpair1.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when(rpair1.forward()).thenReturn(r1f);
        when(rpair1.reverse()).thenReturn(r1r);

        // This is our first read pair so of course it was not seen before
        boolean result = ddmap.hasSeen(rpair1);
        assertFalse(result);

        // Now add a duplicate read pair and make sure it is recognized as a duplicate
        ReadPair rpair2 = Mockito.mock(ReadPair.class);

        when(rpair2.isTrans()).thenReturn(false);
        when(rpair2.getFivePrimeEndPosOfR1()).thenReturn(10);
        when(rpair2.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when(rpair2.forward()).thenReturn(r2f);
        when(rpair2.reverse()).thenReturn(r2r);

        result = ddmap.hasSeen(rpair2);
        assertTrue(result);

        // we have performed two queries
        assertEquals(2, ddmap.getNumOfQueries());
        // There is only one fake chromosome
        assertEquals(1, ddmap.getNumOfChrPairKeys());

    }
        // Now test two read pairs that have the same coordinates (5' end positions) but different orientations that cannot be explained by a duplication event
        @Test
        void recognizeDuplicateTestDifferingOrientation() {
        boolean useRelativeOrientation=true;
        DeDupMap ddmap2 = new DeDupMap(useRelativeOrientation);

        ReadPair rpair3 = Mockito.mock(ReadPair.class);
        SAMRecord r3f=Mockito.mock(SAMRecord.class);
        SAMRecord r3r=Mockito.mock(SAMRecord.class);
        when(rpair3.isTrans()).thenReturn(false);
        when (rpair3.getFivePrimeEndPosOfR1()).thenReturn(10);
        when (rpair3.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when (rpair3.forward()).thenReturn(r3f);
        when (rpair3.reverse()).thenReturn(r3r);
        when (r3f.getReferenceName()).thenReturn("chr1");
        when (r3r.getReferenceName()).thenReturn("chr1");
        when (rpair3.getRelativeOrientationTag()).thenReturn("F1F2");

        // add another read pair that differs with respect to orientation only
        ReadPair rpair4 = Mockito.mock(ReadPair.class);
        SAMRecord r4f=Mockito.mock(SAMRecord.class);
        SAMRecord r4r=Mockito.mock(SAMRecord.class);
        when(rpair4.isTrans()).thenReturn(false);
        when (rpair4.getFivePrimeEndPosOfR1()).thenReturn(10);
        when (rpair4.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when (rpair4.forward()).thenReturn(r4f);
        when (rpair4.reverse()).thenReturn(r4r);
        when (r4f.getReferenceName()).thenReturn("chr1");
        when (r4r.getReferenceName()).thenReturn("chr1");
        when (rpair4.getRelativeOrientationTag()).thenReturn("F1R2"); // two read pairs with F1F2 and F1R2 cannot originate from duplication event

        // query and insert rpair3
        boolean result = ddmap2.hasSeen(rpair3);
        assertFalse(result); // rpair3 should be unknown for DeDupMap

        // query and insert rpair4
        result = ddmap2.hasSeen(rpair4);
        assertFalse(result); // this read pair should be inserted, because of the different orientation
    }

    /**
     * Test compares removal of duplicates with and without consideration of relative orientation of reads.
     * Eight read pairs each with the same five prime end positions but different orientations are created
     * and inserted into the DeDupMap. For the map that does not take into account relative orientation,
     * we expect only one insertion, whereas for the map that does take into account orientation,
     * we expect four insertions, because there are four pairs of read pairs that have equivalent orientation.
     * See here for an explanatory figure: src/test/resources/data/testDeDup/doc/deDupDoc.pdf
     */
    @Test
    void testDuplicateRemovalWithAndWithoutConsiderationOfOrienation() {

        DeDupMap ddmap = new DeDupMap(false);
        DeDupMap ddmap_ori = new DeDupMap(true);


        // first pair of duplicated pairs
        // ------------------------------

        ReadPair rpair1 = Mockito.mock(ReadPair.class);
        when(rpair1.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair1.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair1.getFivePrimeEndPosOfR1()).thenReturn(1000);
        when(rpair1.getFivePrimeEndPosOfR2()).thenReturn(2000);
        when(rpair1.getRelativeOrientationTag()).thenReturn("F1R2");

        ReadPair rpair2 = Mockito.mock(ReadPair.class);
        when(rpair2.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair2.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair2.getFivePrimeEndPosOfR1()).thenReturn(2000);
        when(rpair2.getFivePrimeEndPosOfR2()).thenReturn(1000);
        when(rpair2.getRelativeOrientationTag()).thenReturn("F2R1"); // only difference to rpair1

        ddmap.hasSeen(rpair1);
        ddmap.hasSeen(rpair2);

        ddmap_ori.hasSeen(rpair1);
        ddmap_ori.hasSeen(rpair2);

        assertEquals(1,ddmap.getNumOfInsertions());
        assertEquals(1,ddmap_ori.getNumOfInsertions());


        // second pair of duplicated pairs
        // -------------------------------

        ReadPair rpair3 = Mockito.mock(ReadPair.class);
        when(rpair3.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair3.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair3.getFivePrimeEndPosOfR1()).thenReturn(1000);
        when(rpair3.getFivePrimeEndPosOfR2()).thenReturn(2000);
        when(rpair3.getRelativeOrientationTag()).thenReturn("R1F2");

        ReadPair rpair4 = Mockito.mock(ReadPair.class);
        when(rpair4.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair4.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair4.getFivePrimeEndPosOfR1()).thenReturn(2000);
        when(rpair4.getFivePrimeEndPosOfR2()).thenReturn(1000);
        when(rpair4.getRelativeOrientationTag()).thenReturn("R2F1"); // only difference to rpair3

        ddmap.hasSeen(rpair3);
        ddmap.hasSeen(rpair4);

        ddmap_ori.hasSeen(rpair3);
        ddmap_ori.hasSeen(rpair4);

        assertEquals(1, ddmap.getNumOfInsertions());
        assertEquals(2, ddmap_ori.getNumOfInsertions());


        // third pair of duplicated pairs
        // ------------------------------

        ReadPair rpair5 = Mockito.mock(ReadPair.class);
        when(rpair5.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair5.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair5.getFivePrimeEndPosOfR1()).thenReturn(1000);
        when(rpair5.getFivePrimeEndPosOfR2()).thenReturn(2000);
        when(rpair5.getRelativeOrientationTag()).thenReturn("F1F2");

        ReadPair rpair6 = Mockito.mock(ReadPair.class);
        when(rpair6.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair6.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair6.getFivePrimeEndPosOfR1()).thenReturn(2000);
        when(rpair6.getFivePrimeEndPosOfR2()).thenReturn(1000);
        when(rpair6.getRelativeOrientationTag()).thenReturn("F2F1"); // only difference to rpair5

        ddmap.hasSeen(rpair5);
        ddmap.hasSeen(rpair6);

        ddmap_ori.hasSeen(rpair5);
        ddmap_ori.hasSeen(rpair6);

        assertEquals(1, ddmap.getNumOfInsertions());
        assertEquals(3, ddmap_ori.getNumOfInsertions());

        // fourth pair of duplicated pairs
        // -------------------------------

        ReadPair rpair7 = Mockito.mock(ReadPair.class);
        when(rpair7.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair7.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair7.getFivePrimeEndPosOfR1()).thenReturn(1000);
        when(rpair7.getFivePrimeEndPosOfR2()).thenReturn(2000);
        when(rpair7.getRelativeOrientationTag()).thenReturn("R1R2");

        ReadPair rpair8 = Mockito.mock(ReadPair.class);
        when(rpair8.getReferenceSequenceOfR1()).thenReturn("chr1");
        when(rpair8.getReferenceSequenceOfR2()).thenReturn("chr1");
        when(rpair8.getFivePrimeEndPosOfR1()).thenReturn(2000);
        when(rpair8.getFivePrimeEndPosOfR2()).thenReturn(1000);
        when(rpair8.getRelativeOrientationTag()).thenReturn("R2R1"); // only difference to rpair7

        ddmap.hasSeen(rpair7);
        ddmap.hasSeen(rpair8);

        ddmap_ori.hasSeen(rpair7);
        ddmap_ori.hasSeen(rpair8);

        assertEquals(1, ddmap.getNumOfInsertions());
        assertEquals(4, ddmap_ori.getNumOfInsertions());
    }
}
