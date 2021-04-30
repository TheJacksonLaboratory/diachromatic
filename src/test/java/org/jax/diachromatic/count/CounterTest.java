package org.jax.diachromatic.count;

import htsjdk.samtools.SamReader;
import org.jax.diachromatic.align.Digest;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.align.DigestPair;
import org.jax.diachromatic.align.ReadPair;
import org.junit.jupiter.api.Test;
import static org.junit.Assert.assertEquals;
import org.mockito.Mockito;
import static org.mockito.Mockito.*;

public class CounterTest {

    /**
     * Test incrementing of interaction counts.
     */
    @Test
    public void testIncrementDigestPair() {

        Digest d1 = Mockito.mock(Digest.class);
        Digest d2 = Mockito.mock(Digest.class);

        when(d1.getDigestStartPosition()).thenReturn(1000);
        when(d2.getDigestStartPosition()).thenReturn(2000);

        DigestPair dp1 = Mockito.mock(DigestPair.class);

        when(dp1.forward()).thenReturn(d1);
        when(dp1.reverse()).thenReturn(d2);

        ReadPair rp1 = Mockito.mock(ReadPair.class);
        when(rp1.isTwisted()).thenReturn(false);

        ReadPair rp2 = Mockito.mock(ReadPair.class);
        when(rp2.isTwisted()).thenReturn(true);

        DigestMap dm = Mockito.mock(DigestMap.class);
        SamReader samReader = Mockito.mock(SamReader.class);
        Counter counter = new Counter(samReader, dm, "fakeOutputPathWithPrefix",false);

        // test whether the same count is incremented for different orders of digests in a given pair.
        counter.incrementDigestPair(dp1,rp1);
        when(dp1.forward()).thenReturn(d2);
        when(dp1.reverse()).thenReturn(d1);
        counter.incrementDigestPair(dp1,rp1);
        assertEquals(1,counter.getInteractionCount());

        // test whether incrementing with new digest pair object consisting of the same digests increases the interaction count
        // TODO: It does! Discuss with Nick if this is intended behavior.
        DigestPair dp2 = Mockito.mock(DigestPair.class);
        when(dp2.forward()).thenReturn(d1);
        when(dp2.reverse()).thenReturn(d2);
        counter.incrementDigestPair(dp2,rp1);
        assertEquals(2,counter.getInteractionCount());

        // test whether simple and twisted read pairs are counted correctly
        counter.incrementDigestPair(dp1,rp2);
        assertEquals(2,counter.getSimpleTwistedCountForDigestPair(dp1).simple_1 + counter.getSimpleTwistedCountForDigestPair(dp1).simple_2);
        assertEquals(1,counter.getSimpleTwistedCountForDigestPair(dp1).twisted_1 + counter.getSimpleTwistedCountForDigestPair(dp1).twisted_2);
    }
}