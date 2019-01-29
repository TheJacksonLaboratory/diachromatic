package org.jax.diachromatic.allelespec;

import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

/**
 * This class stores information about loci that are not
 * homozygous wildtype or alt, and therefore can be regarded
 * as being candidates for allele-specific loci.
 * The implementation for now is pretty naive -- TODO explore how to be more efficient
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 */
public class AlleleSpecificLocus {

    private final String chromosome;
    private final int position;
    private final int total;
    private final int count1;
    private final int count2;



    public AlleleSpecificLocus(String chrom, int pos, int total, int count1, int count2){
        this.chromosome=chrom;
        this.position=pos;
        this.total=total;
        this.count1=count1;
        this.count2=count2;
    }


    public double binomialTestPvalue() {
        double probability=0.5000;
        return new BinomialTest().binomialTest(total, count2,  probability, AlternativeHypothesis.TWO_SIDED);
    }




    public String getChromosome() {
        return chromosome;
    }

    public int getPosition() {
        return position;
    }

    public int getTotal() {
        return total;
    }

    public int getCount1() {
        return count1;
    }

    public int getCount2() {
        return count2;
    }

    @Override
    public String toString() {
        return String.format("%s:%s count 1: %d count 2: %d",chromosome,position,count1,count2);
    }
}
