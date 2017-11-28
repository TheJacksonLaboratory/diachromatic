package org.jax.diachromatic.map;

import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class SamBitflagFilterTest {

    /**
     * SAM flag 4 means read unmapped (0x4)
     * SAM flag 0 means no flag set, read must be mapped
     */
    @Test
    public void testUnmappedSingleEnd() {
        int flag=4;
        assertTrue(SamBitflagFilter.segmentIsUnmapped(flag));
        flag = 0;
        assertFalse(SamBitflagFilter.segmentIsUnmapped(flag));
    }



}
