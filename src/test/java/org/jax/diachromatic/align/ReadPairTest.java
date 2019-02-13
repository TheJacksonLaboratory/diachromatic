package org.jax.diachromatic.align;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMRecord;
import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.Test;
import org.mockito.Mockito;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.mockito.Mockito.when;

public class ReadPairTest {


    @Test
    public void testMultiMappedReads() throws DiachromaticException {
        SAMRecord record1 = Mockito.mock(SAMRecord.class);
        SAMRecord record2 = Mockito.mock(SAMRecord.class);
        Map<String, List<Digest>> emptymap=ImmutableMap.of();
        DigestMap mockDigestMap = Mockito.mock(DigestMap.class);
        int upperFragSize = 250;
        int lowerFragSize = 200;
        boolean stringent = false;

        when(record1.getReferenceName()).thenReturn("fakeChromosome1");
        when(record2.getReferenceName()).thenReturn("fakeChromosome2");

        when(record1.getReadUnmappedFlag()).thenReturn(false);
        when(record2.getReadUnmappedFlag()).thenReturn(false);

        // test 1
        when(record1.getMappingQuality()).thenReturn(31); // mapping quality is sufficient
        when(record1.getAttribute("AS")).thenReturn(0);
        when(record1.getAttribute("XS")).thenReturn(-5); // but the distance between the scores of the best and second best alignment is two small
        ReadPair rpair = new ReadPair(record1,record2,mockDigestMap,stringent);
        assertFalse(rpair.isPaired());

        // test 2
        when(record1.getAttribute("AS")).thenReturn(0);
        when(record1.getAttribute("XS")).thenReturn(-20); // distance between the scores of the best and second best alignment is sufficiently large
        when(record1.getMappingQuality()).thenReturn(25); // but the mapping quality is insufficient
        rpair = new ReadPair(record1,record2,mockDigestMap,stringent);
        assertFalse(rpair.isPaired());

        // test 3 now use the stringent definition of multi mapping
        stringent = true;
        when(record1.getMappingQuality()).thenReturn(31); // mapping quality is sufficient
        when(record1.getAttribute("AS")).thenReturn(0);
        when(record1.getAttribute("XS")).thenReturn(-20); // distance between the scores of the best and second best alignment is sufficiently large
        rpair = new ReadPair(record1,record2,mockDigestMap,stringent);
        assertFalse(rpair.isPaired()); // but this is not enough for the stringent mode; the reads remain unpaired

        // positive cannot be tested using incomplete fake mock objects
    }

    @Test
    public void testUniMappedReads() throws DiachromaticException {
        SAMRecord record1 = Mockito.mock(SAMRecord.class);
        SAMRecord record2 = Mockito.mock(SAMRecord.class);
        Map<String, List<Digest>> emptymap=ImmutableMap.of();
        DigestMap mockDigestMap = Mockito.mock(DigestMap.class);

        when(record1.getAlignmentStart()).thenReturn(100);
        when(record1.getAlignmentEnd()).thenReturn(200);

        when(record2.getAlignmentStart()).thenReturn(400);
        when(record2.getAlignmentEnd()).thenReturn(500);


        Digest digest = Mockito.mock(Digest.class);
        DigestPair mockDigestPair = Mockito.mock(DigestPair.class);
        when(mockDigestPair.forward()).thenReturn(digest);
        when(mockDigestPair.reverse()).thenReturn(digest);

        when(mockDigestMap.getDigestPair("fakeChromosome1", 200, "fakeChromosome2", 400)).thenReturn(mockDigestPair);
        int upperFragSize = 250;
        int lowerFragSize = 200;
        boolean stringent = false;

        when(record1.getReferenceName()).thenReturn("fakeChromosome1");
        when(record2.getReferenceName()).thenReturn("fakeChromosome2");

        when(record1.getReadUnmappedFlag()).thenReturn(false);
        when(record2.getReadUnmappedFlag()).thenReturn(false);

        when(record1.getAttribute("XS")).thenReturn(null);
        when(record2.getAttribute("XS")).thenReturn(null);

        when(record1.getReadNegativeStrandFlag()).thenReturn(true);
        when(record2.getReadNegativeStrandFlag()).thenReturn(false);

        ReadPair rpair = new ReadPair(record1,record2,mockDigestMap,stringent);
        assertTrue(rpair.isPaired());
    }
}
