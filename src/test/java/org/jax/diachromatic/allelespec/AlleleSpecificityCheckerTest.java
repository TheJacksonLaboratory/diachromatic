package org.jax.diachromatic.allelespec;

import org.junit.Test;

public class AlleleSpecificityCheckerTest{
    /** Just to get started, a local file. */
    private static String bamPath="/Users/peterrobinson/Documents/data/diachromatic/HG00158_chr1_chr22_toy_example.bam";


    @Test
    public void test1() {
        AlleleSpecificityChecker checker = new AlleleSpecificityChecker(bamPath,"");
        checker.readAlleles();
    }


}
