package org.jax.diachromatic.allelespec;

import org.junit.BeforeClass;
import org.junit.Test;
import org.junit.jupiter.api.Disabled;

import java.io.FileNotFoundException;
import java.net.URL;

import static org.junit.Assert.assertEquals;
@Disabled("Only used to output text during development.")
public class AlleleSpecificityCheckerTest{
    // Path to BAM file created with help of {@link MakeFilesForTest}.
    private static String bamPath;

    @BeforeClass
    public static void init() throws FileNotFoundException {
        ClassLoader classLoader = AlleleSpecificityCheckerTest.class.getClassLoader();
        URL resource = classLoader.getResource("data/allelespec/sample.sorted.bam");
        if (resource==null){
            throw new FileNotFoundException("Could not find sample.sorted.bam file");
        }
        bamPath = resource.getFile();
    }

    @Test
    public void test1() {
        AlleleSpecificityChecker checker = new AlleleSpecificityChecker(bamPath,"");
        checker.readAlleles();
        assertEquals(2,checker.getLocuslist().size());
//        for (AlleleSpecificLocus loc : checker.getLocuslist()) {
//            System.err.println(loc);
//            System.err.println("p="+loc.binomialTestPvalue());
//        }
    }
}
