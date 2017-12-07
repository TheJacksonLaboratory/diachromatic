package org.jax.diachromatic.map;

/**
 * This class provides a set of static functions to perform SAM format bitflag filtering
 * for reads that did not map or had other problems.
 */
public class SamBitflagFilter {
    /* multiple segments, e.g., a read inputSAMfiles. The bit is set if the read is part of a paired end inputSAMfiles. */
    public static final int TEMPLATE_HAS_MULTIPLE_SEGMENTS = 0x1;

    public static final int EACH_SEGMENT_PROPERLY_ALIGNED = 0x2;

    public static final int SEGMENT_IS_UNMAPPED = 0x4;

    public static final int MATE_IS_UNMAPPED = 0x8;

    public static final int SEG_IS_REVCOMP= 0x10;

    public static final int SEQ_NEXT_SEGMENT_IS_REVCOMP = 0x20;

    public static final int FIRST_SEGMENT_IN_TEMPLATE = 0x40;

    public static final int LAST_SEGMENT_IN_TEMPLATE = 0x80;

    public static final int SECONDARY_ALIGNMENT=0x100;

    public static final int SEGMENT_DOES_NOT_PASS_QC = 0x200;

    public static final int SEGMENT_IS_PCR_OPTICAL_DUPLICATE = 0x400;

    public static final int SUPPLEMENTARY_ALIGNMENT = 0x800;


//    public static boolean segmentIsUnmapped(int flag) {
//        return (flag & SEGMENT_IS_UNMAPPED)== SEGMENT_IS_UNMAPPED;
//    }
//


    public static  void debugDisplayBitflag(int flag) {
        System.out.println("flag: "+flag);
        if ((flag & TEMPLATE_HAS_MULTIPLE_SEGMENTS)==TEMPLATE_HAS_MULTIPLE_SEGMENTS) {
            System.out.println("\t0x1; template has multiple segments");
        }
        if ((flag & EACH_SEGMENT_PROPERLY_ALIGNED)==EACH_SEGMENT_PROPERLY_ALIGNED) {
            System.out.println("\t0x2; EACH_SEGMENT_PROPERLY_ALIGNED");
        }
        if ((flag & SEGMENT_IS_UNMAPPED)==SEGMENT_IS_UNMAPPED) {
            System.out.println("\t0x4; SEGMENT_IS_UNMAPPED");
        }
        if ((flag & MATE_IS_UNMAPPED)==MATE_IS_UNMAPPED) {
            System.out.println("\t0x8; MATE_IS_UNMAPPED");
        }
        if ((flag & SEG_IS_REVCOMP)==SEG_IS_REVCOMP) {
            System.out.println("\t0x10; SEG_IS_REVCOMP");
        }
        if ((flag & SEQ_NEXT_SEGMENT_IS_REVCOMP)==SEQ_NEXT_SEGMENT_IS_REVCOMP) {
            System.out.println("\t0x20; SEQ_NEXT_SEGMENT_IS_REVCOMP");
        }
        if ((flag & FIRST_SEGMENT_IN_TEMPLATE)==FIRST_SEGMENT_IN_TEMPLATE) {
            System.out.println("\t0x40; FIRST_SEGMENT_IN_TEMPLATE");
        }
        if ((flag & LAST_SEGMENT_IN_TEMPLATE)==LAST_SEGMENT_IN_TEMPLATE) {
            System.out.println("\t0x80; LAST_SEGMENT_IN_TEMPLATE");
        }
        if ((flag & SECONDARY_ALIGNMENT)==SECONDARY_ALIGNMENT) {
            System.out.println("\t0x100; SECONDARY_ALIGNMENT");
        }
        if ((flag & SEGMENT_DOES_NOT_PASS_QC)==SEGMENT_DOES_NOT_PASS_QC) {
            System.out.println("\t0x200; SEGMENT_DOES_NOT_PASS_QC");
        }
        if ((flag & SEGMENT_IS_PCR_OPTICAL_DUPLICATE)==SEGMENT_IS_PCR_OPTICAL_DUPLICATE) {
            System.out.println("\t0x400; SEGMENT_IS_PCR_OPTICAL_DUPLICATE");
        }
        if ((flag & SUPPLEMENTARY_ALIGNMENT)==SUPPLEMENTARY_ALIGNMENT) {
            System.out.println("\t0x800; SUPPLEMENTARY_ALIGNMENT");
        }


    }




}
