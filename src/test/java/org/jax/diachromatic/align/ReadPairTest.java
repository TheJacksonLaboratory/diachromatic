package org.jax.diachromatic.align;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.SAMRecord;
import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.Test;
import org.mockito.Mockito;

import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertFalse;
import static org.mockito.Mockito.when;

public class ReadPairTest {


    @Test
    public void test1() throws DiachromaticException {
        SAMRecord record1 = Mockito.mock(SAMRecord.class);
        SAMRecord record2 = Mockito.mock(SAMRecord.class);
        Map<String, List<Digest>> emptymap=ImmutableMap.of();
        DigestMap mockDigestMap = Mockito.mock(DigestMap.class);
        int upperFrage = 250;
        int lowerfrag=200;
        boolean stringent = false;

        when(record1.getReadUnmappedFlag()).thenReturn(false);
        when(record2.getReadUnmappedFlag()).thenReturn(false);
        when(record1.getAttribute("XS")).thenReturn(0);
        when(record1.getAttribute("AS")).thenReturn(-5);
        when(record2.getAttribute("XS")).thenReturn(null);
        when(record1.getReferenceName()).thenReturn("fakeChromosome1");
        when(record2.getReferenceName()).thenReturn("fakeChromosome2");

        when(record1.getMappingQuality()).thenReturn(20);

        ReadPair rpair = new ReadPair(record1,record2,emptymap,mockDigestMap,upperFrage,lowerfrag,stringent);


        assertFalse(rpair.isPaired());

        // now test negative
        when(record1.getAttribute("XS")).thenReturn(-15);
        when(record1.getAttribute("AS")).thenReturn(-20);

        rpair = new ReadPair(record1,record2,emptymap,mockDigestMap,upperFrage,lowerfrag,stringent);
        assertFalse(rpair.isPaired());

    }



}
