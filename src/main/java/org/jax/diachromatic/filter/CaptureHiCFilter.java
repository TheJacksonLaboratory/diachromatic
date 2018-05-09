package org.jax.diachromatic.filter;

/**
 * (hicup: Since valid Hi-C reads do not represent one continuous genomic sequence, the pipeline maps forward and
 * reverse reads independently and then re-pairs sequences in which both ends aligned unambiguously to the genome.)
 * TODO Do we want to align from within the Java process with bowtie2 ?
 */
public class CaptureHiCFilter {

    public CaptureHiCFilter() {}
}
