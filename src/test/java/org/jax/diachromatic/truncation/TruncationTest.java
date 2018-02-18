package org.jax.diachromatic.truncation;

import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.junit.BeforeClass;
import org.junit.ClassRule;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

/**
 * A test dataset was prepared. There is a suffix on the sequence name that indicates the expected length of the
 * fragment following the truncation step. For instance
 * @SRR071233.15278601.35_bp NRTG514-16_0001:3:117:15170:5880 length=40
 * For this sequence, the original length is 40nt, and the length after truncation is 35.
 * We perform the truncation step and write the results to a temp file; we then check the length of the
 * truncated sequences.
 */
public class TruncationTest {

    @ClassRule
    public static TemporaryFolder tempFolder= new TemporaryFolder();
    /** Key -- sequence name line value--the sequence (for FASTQ 1) */
    private static Map<String,String> fastq1map=new HashMap<>();
    /** Key -- sequence name line value--the sequence (for FASTQ 2) */
    private static Map<String,String> fastq2map=new HashMap<>();




    @BeforeClass
    public static void init() throws IOException {
        ClassLoader classLoader = TruncatorTest.class.getClassLoader();
        String fastq_1 = classLoader.getResource("data/fastq/hg19_HindIII_test_data_truncation_1.fastq").getFile();
        String fastq_2 = classLoader.getResource("data/fastq/hg19_HindIII_test_data_truncation_2.fastq").getFile();
        RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
        List<RestrictionEnzyme> enzymelist = new ArrayList<>();
        enzymelist.add(hindIII);
               File tmpFolder = tempFolder.newFolder();
        Truncator truncator= new Truncator(tmpFolder.getAbsolutePath(),fastq_1,fastq_2,hindIII,"test");
        System.err.println(tmpFolder.getAbsolutePath());
        truncator.parseFASTQ();
        String fastQ1=truncator.getOutputFastqPath_1();
        String fastQ2=truncator.getOutputFastqPath_2();
        BufferedReader br = new BufferedReader(new FileReader(fastQ1));
        String line;
        int i=0;
        String currentName=null;
        while ((line=br.readLine())!=null) {
            if (i%4==0) {
                // name line
                int j=line.indexOf("_bp")+3;
                currentName=line.substring(0,j);
            } else if (i%4==1) {
                // sequence line
                fastq1map.put(currentName,line);
            }
            i++;
        }
        br.close();
        br = new BufferedReader(new FileReader(fastQ2));
        i=0;
        currentName=null;
        while ((line=br.readLine())!=null) {
            if (i%4==0) {
                // name line
                int j=line.indexOf(" ");
                currentName=line.substring(0,j);
            } else if (i%4==1) {
                // sequence line
                fastq2map.put(currentName,line);
            }
            i++;
        }

    }

    /** The first read pair. Only first read is truncated to 35 bp */
    @Test
    public void testSRR071233_15278601_35_bp() {
        String seq1=fastq1map.get("@SRR071233.15278601.35_bp");
        assertNotNull(seq1);
        int expected=35;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.15278601.40_bp");
        assertNotNull(seq2);
        expected=40;
        assertEquals(expected,seq2.length());
    }

    /** The second read pair. Only second read is truncated to 29 bp */
    @Test public void testSRR071233_15562685_40_bp () {
        String seq1=fastq1map.get("@SRR071233.15562685.40_bp");
        assertNotNull(seq1);
        int expected=40;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.15562685.29_bp");
        assertNotNull(seq2);
        expected=29;
        assertEquals(expected,seq2.length());
    }

    /** The third read pair. Only first read is truncated to 32 bp */
    @Test public void testSRR071233_12020325_32_bp() {
        String seq1 = fastq1map.get("@SRR071233.12020325.32_bp");
        assertNotNull(seq1);
        int expected=32;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.12020325.40_bp");
        assertNotNull(seq2);
        expected=40;
        assertEquals(expected,seq2.length());
    }

    /** The fourth read pair. Only second read is truncated to 24 bp */
    @Test public void testSRR071233_13947609_40_bp() {
        String seq1 = fastq1map.get("@SRR071233.13947609.40_bp");
        assertNotNull(seq1);
        int expected=40;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.13947609.24_bp");
        assertNotNull(seq2);
        expected=24;
        assertEquals(expected,seq2.length());
    }

    /** The fifth read pair. Only first read is truncated to 25 bp */
    @Test public void testSRR071233_14697319_25_bp() {
        String seq1 = fastq1map.get("@SRR071233.14697319.25_bp");
        assertNotNull(seq1);
        int expected=25;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.14697319.40_bp");
        expected=40;
        assertEquals(expected,seq2.length());
    }
    /** The sixth read pair. Only second read is truncated to 27 bp */
    @Test public void testSRR071233_12307569_40_bp() {
        String seq1 = fastq1map.get("@SRR071233.12307569.40_bp");
        assertNotNull(seq1);
        int expected=40;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.12307569.27_bp");
        expected=27;
        assertEquals(expected,seq2.length());
    }
    /** The seventh read pair. Only first read is truncated to 27 bp */
    @Test public void testSRR071233_15497273_27_bp() {
        String seq1 = fastq1map.get("@SRR071233.15497273.27_bp");
        assertNotNull(seq1);
        int expected=27;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.15497273.40_bp");
        expected=40;
        assertEquals(expected,seq2.length());
    }

    /** The eigth read pair. Only second read is truncated to 23 bp */
    @Test public void testSRR071233_13860235_40_bp() {
        String seq1 = fastq1map.get("@SRR071233.13860235.40_bp");
        assertNotNull(seq1);
        int expected=40;
        assertEquals(expected,seq1.length());
        String seq2=fastq2map.get("@SRR071233.13860235.23_bp");
        expected=23;
        assertEquals(expected,seq2.length());
    }

}
