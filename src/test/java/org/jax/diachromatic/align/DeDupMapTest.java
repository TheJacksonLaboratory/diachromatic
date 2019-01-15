package org.jax.diachromatic.align;


import htsjdk.samtools.SAMRecord;
import org.junit.jupiter.api.Test;
import org.mockito.Mockito;


import static org.junit.jupiter.api.Assertions.*;
import static org.mockito.Mockito.when;

class DeDupMapTest {


    @Test
    void recognizeDuplicateTest() {
        boolean useRelativeOrientation=false;
        DeDupMap ddmap = new DeDupMap(useRelativeOrientation);
        ReadPair rpair1 = Mockito.mock(ReadPair.class);
        // mock SAMRecords, forward and reverse for the two read pairs.
        SAMRecord r1f=Mockito.mock(SAMRecord.class);
        SAMRecord r1r=Mockito.mock(SAMRecord.class);
        SAMRecord r2f=Mockito.mock(SAMRecord.class);
        SAMRecord r2r=Mockito.mock(SAMRecord.class);
        when (r1f.getReferenceName()).thenReturn("chrZ");
        when (r1r.getReferenceName()).thenReturn("chrZ");
        when (r2f.getReferenceName()).thenReturn("chrZ");
        when (r2r.getReferenceName()).thenReturn("chrZ");

        when(rpair1.isTrans()).thenReturn(false);
        when (rpair1.getFivePrimeEndPosOfR1()).thenReturn(10);
        when (rpair1.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when (rpair1.forward()).thenReturn(r1f);
        when (rpair1.reverse()).thenReturn(r1r);

        // This is our first readpair so of course it was not seen before.
        boolean result = ddmap.hasSeen(rpair1);
        assertFalse(result);

        // Now add a duplicate readpair and make sure it is recognized as a duplicate
        ReadPair rpair2 = Mockito.mock(ReadPair.class);

        when(rpair2.isTrans()).thenReturn(false);
        when (rpair2.getFivePrimeEndPosOfR1()).thenReturn(10);
        when (rpair2.getFivePrimeEndPosOfR2()).thenReturn(1010);
        when (rpair2.forward()).thenReturn(r2f);
        when (rpair2.reverse()).thenReturn(r2r);

        result = ddmap.hasSeen(rpair2);
        assertTrue(result);

        // we have performed two queries
        assertEquals(2,ddmap.getNumOfQueries());
        // There is only one fake chromosome
        assertEquals(1,ddmap.getNumOfChrPairKeys());


    }



}
