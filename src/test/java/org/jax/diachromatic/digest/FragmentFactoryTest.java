package org.jax.diachromatic.digest;

import org.jax.diachromatic.exception.DiachromaticException;
import org.junit.BeforeClass;
import org.junit.ClassRule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class FragmentFactoryTest {

    @ClassRule
    public static TemporaryFolder folder= new TemporaryFolder();

    private static String pathToTempFile=null;

    private static String genomeDirectoryPath=null;
    private static FragmentFactory factory=null;

    @BeforeClass
    public static void setup() throws IOException {
        ClassLoader classLoader = FragmentFactoryTest.class.getClassLoader();
        String fastaPath = classLoader.getResource("data/chrUn_KI270745v1.fa").getFile();
        genomeDirectoryPath = classLoader.getResource("data").getFile();
        File tempFile = folder.newFile("file.txt");
        pathToTempFile=tempFile.getAbsolutePath();
        factory=new FragmentFactory(genomeDirectoryPath, tempFile.getAbsolutePath());
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


    /** We will test cutting our test file chrUn_KI270745v1.fa
     * with HindIII. We used the hiCUP digester as follows
     * <pre>
     *     $ ./hicup_digester --re1 A^AGCTT -genome hg38 chrUn_KI270745v1.fa
     * </pre>
     * Note: hindIII = new RestrictionEnzyme("HindIII", "A^AGCTT");
     * This produced
     * <pre>
     * Genome:hg38	Restriction_Enzyme1:re1_unspecified [A^AGCTT]	Restriction_Enzyme2:None	Hicup digester version 0.5.10
     Chromosome	Fragment_Start_Position	Fragment_End_Position	Fragment_Number	RE1_Fragment_Number	5'_Restriction_Site	3'_Restriction_Site
     chrUn_KI270745v1	1	8241	1	1	None	Re1
     chrUn_KI270745v1	8242	9365	2	2	Re1	Re1
     chrUn_KI270745v1	9366	12175	3	3	Re1	Re1
     chrUn_KI270745v1	12176	12399	4	4	Re1	Re1
     chrUn_KI270745v1	12400	15584	5	5	Re1	Re1
     chrUn_KI270745v1	15585	15808	6	6	Re1	Re1
     chrUn_KI270745v1	15809	19006	7	7	Re1	Re1
     chrUn_KI270745v1	19007	19229	8	8	Re1	Re1
     chrUn_KI270745v1	19230	29179	9	9	Re1	Re1
     chrUn_KI270745v1	29180	29246	10	10	Re1	Re1
     chrUn_KI270745v1	29247	29818	11	11	Re1	Re1
     chrUn_KI270745v1	29819	35094	12	12	Re1	Re1
     chrUn_KI270745v1	35095	38321	13	13	Re1	Re1
     chrUn_KI270745v1	38322	38683	14	14	Re1	Re1
     chrUn_KI270745v1	38684	41891	15	15	Re1	None
     * </pre>
     */
    @Test
    public void testHindIIIcut() throws IOException {
        System.err.println("Indexing FASTA files");
        // 1. cut with DpnII
        try {
            factory.indexFASTAfilesIfNeeded();
            String testEnzyme = "HindIII";
            List<String> enzymes=new ArrayList<>();
            enzymes.add(testEnzyme);
            factory.digestGenome(enzymes);
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }
        // input the temporary file and check the results
        BufferedReader br = new BufferedReader(new FileReader(pathToTempFile));
        String line=null;
        List<String> lineList=new ArrayList<>();
        while ((line=br.readLine())!=null) {
            lineList.add(line);
        }
        br.close();
        // Go through the lines and test equality to hicup output
        // Note we have slightly modified the format
        line = lineList.get(1);
        String A[];
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","1","8241","1","None","HindIII"));
        line = lineList.get(2);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","8242","9365","2","HindIII","HindIII"));
        line = lineList.get(3);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","9366","12175","3","HindIII","HindIII"));
        line = lineList.get(4);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","12176","12399","4","HindIII","HindIII"));
        line = lineList.get(5);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","12400","15584","5","HindIII","HindIII"));
        line = lineList.get(6);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","15585","15808","6","HindIII","HindIII"));
        line = lineList.get(7);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","15809","19006","7","HindIII","HindIII"));
        line = lineList.get(8);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","19007","19229","8","HindIII","HindIII"));
        line = lineList.get(9);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","19230","29179","9","HindIII","HindIII"));
        line = lineList.get(10);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","29180","29246","10","HindIII","HindIII"));
        line = lineList.get(11);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","29247","29818","11","HindIII","HindIII"));
        line = lineList.get(12);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","29819","35094","12","HindIII","HindIII"));
        line = lineList.get(13);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","35095","38321","13","HindIII","HindIII"));
        line = lineList.get(14);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","38322","38683","14","HindIII","HindIII"));
        line = lineList.get(15);
        A=line.split("\t");
        assertTrue(testline(A,"chrUn_KI270745v1","38684","41891","15","HindIII","None"));
    }

    /**
     * This function is intended to be used together with {@link #testHindIIIcut()}, and it tests whether
     * a line of the output file of digestion is as expected
     * @param A String array with 6 fields
     * @param chrom chromosome
     * @param start start position of restriction fragment
     * @param end end position of restriction fragment
     * @param nr number of restriction fragment
     * @param enz1 restriction enzyme on 5' side
     * @param enz2 restriction enzyme on 3' side
     * @return true if data in A matches the other parameters, otherwise false.
     */
    private boolean testline(String[] A, String chrom, String start, String end, String nr, String enz1, String enz2) {
        if (A.length!=6) return false;
        return (A[0].equals(chrom) &&
        A[1].equals(start) &&
        A[2].equals(end) &&
        A[3].equals(nr) &&
        A[4].equals(enz1) &&
        A[5].equals(enz2));
    }


}
