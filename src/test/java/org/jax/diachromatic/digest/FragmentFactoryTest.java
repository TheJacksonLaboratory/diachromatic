package org.jax.diachromatic.digest;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.BeforeClass;
import org.junit.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;

public class FragmentFactoryTest {

    private static String genomeDirectoryPath=null;
    private static FragmentFactory factory=null;

    @BeforeClass
    public static void setup() throws IOException {
        ClassLoader classLoader = FragmentFactoryTest.class.getClassLoader();
        String fastaPath = classLoader.getResource("data/chrUn_KI270745v1.fa").getFile();
        genomeDirectoryPath = classLoader.getResource("data").getFile();

        factory=new FragmentFactory(genomeDirectoryPath);
    }


    @Test
    public void testGetGenomeDirectoryPath() {
        assertEquals(genomeDirectoryPath,factory.getGenomeDirectoryPath());
    }

    /** There is exactly one FASTA file in src/test/resources/data */
    @Test
    public void testGetGenomeFASTAFileCount() {
        int expected=1;
        assertEquals(expected,factory.getGenomeFileCount());
    }


    @Test
    public void testDpnIIcut() {
        System.err.println("Indexing FASTA files");
        factory.indexFASTAfilesIfNeeded();
        String testEnzyme = "DpnII";
        List<String> enzymes=new ArrayList<>();
        enzymes.add(testEnzyme);
        try {
            factory.digestGenome(enzymes);
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }

        // 1. cut with DpnII
        // 2. make a list of Fragment objects
        // 3. compare with out put of hicup



    }



}
