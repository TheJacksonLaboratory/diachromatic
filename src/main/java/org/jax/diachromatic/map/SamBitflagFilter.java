package org.jax.diachromatic.map;

/**
 * This class provides a set of static functions to perform SAM format bitflag filtering
 * for reads that did not map or had other problems.
 */
public class SamBitflagFilter {
    /* multiple segments, e.g., a read pair. The bit is set if the read is part of a paired end pair. */
    public static final int TEMPLATE_HAS_MULTIPLE_SEGMENTS = 0x1;

    public static final int EACH_SEGMENT_PROPERLY_ALIGNED = 0x2;

    public static final int SEGMENT_IS_UNMAPPED = 0x4;

    public static final int MATE_IS_UNMAPPED = 0x8;

    public static final int SEG_IS_REVCOMP= 0x10;

    public static final int SEQ_NEXT_SEGMENT_IS_REVCOMP = 0x20;

    public static final int FIRST_SEGMENT_IN_TEMPLTE = 0x40;

    public static final int LAST_SEGMENT_IN_TEMPLTE = 0x80;

    public static final int SECONDARY_ALIGNMENT=0x100;

    public static final int SEGMENT_DOES_NOT_PASS_QC = 0x200;

    public static final int SEGMENT_IS_PCR_OPTICAL_DUPLICATE = 0x400;

    public static final int SUPPLEMENTARY_ALIGNMENT = 0x800;


    public static boolean segmentIsUnmapped(int flag) {
        return (flag & SEGMENT_IS_UNMAPPED)== SEGMENT_IS_UNMAPPED;
    }
}
