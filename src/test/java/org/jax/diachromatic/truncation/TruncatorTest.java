package org.jax.diachromatic.truncation;

import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class TruncatorTest {

    private static FastqPairParser parser=null;

  @BeforeClass
    public  static void init() {
      /*  #name   site
        BmtI	GCTAG^C
        KpnI	GGTAC^C
        NsiI	ATGCA^T
        PstI	CTGCA^G
        PvuI	CGAT^CG
        SacI	GAGCT^C
        SacII	CCGCG^G
        SphI	GCATG^C
        FaeI	CATG^
        TaiI	ACGT^  */

      ClassLoader classLoader = TruncatorTest.class.getClassLoader();
      String fastq_1 = classLoader.getResource("data/fastq/test1.fastq").getFile();
      String fastq_2 = classLoader.getResource("data/fastq/test2.fastq").getFile();
      RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
      String ligationSequence=Truncator.fillEnd(hindIII);
      parser=new FastqPairParser(fastq_1,fastq_2,ligationSequence);

    }


    @Test
    public void testBglII() {
        //BglII (5'-A^GATCT-3') cuts between A and G. So, the ligation
        //sequence will be: A + GATC + GATC + T = AGATCGATCT.
        RestrictionEnzyme bglII = new RestrictionEnzyme("BglII","A^GATCT");
        String ligationSequence = Truncator.fillEnd(bglII);
        String expected = "AGATCGATCT";
        assertEquals(expected,ligationSequence);
    }

    @Test
    public void testHindIII() {
        //HindIII  (5'-A^AGCTT-3') cuts between A and A. So, the ligation
        //sequence will be: A + AGCT + AGCT + T = AAGCTAGCTT.
        RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
        String ligationSequence = Truncator.fillEnd(hindIII);
        String expected = "AAGCTAGCTT";
        assertEquals(expected,ligationSequence);
    }

    @Test
    public void testDpnII() {
        //DpnII (5'-^GATC-3') cuts before GATC. So, the ligation
        //sequence will be: GATC + GATC = GATCGATC.
        RestrictionEnzyme dpnII = new RestrictionEnzyme("DpnII","^GATC");
        String ligationSequence = Truncator.fillEnd(dpnII);
        String expected = "GATCGATC";
        assertEquals(expected,ligationSequence);
    }

    @Test
    public void testNlaIII() {
        //NlaIII (5'-CATG^-3') cuts after CATG^. So, the ligation
        //sequence will be: CATG + CATG = CATGCATG.
        RestrictionEnzyme nlaIII = new RestrictionEnzyme("NlaIII","CATG^");
        String ligationSequence = Truncator.fillEnd(nlaIII);
        String expected = "CATGCATG";
        assertEquals(expected,ligationSequence);
    }


    @Test
    public void testAatII() {
        //   GACGT^C --ths ligation sequence will be
        // G-ACGT-ACGT-C = GACGTACGTC
        RestrictionEnzyme aatII = new RestrictionEnzyme("AatII","GACGT^C" );
        String ligationSequence = Truncator.fillEnd(aatII);
        String expected = "GACGTACGTC";
        assertEquals(expected,ligationSequence);

    }

    @Test
    public void testApaI () {
        // ApaI	GGGCC^C
        // G-GGCC-GGCC-C = GGGCCGGCCC
        RestrictionEnzyme apaI = new RestrictionEnzyme("ApaI","GGGCC^C" );
        String ligationSequence = Truncator.fillEnd(apaI);
        String expected = "GGGCCGGCCC";
        assertEquals(expected,ligationSequence);
    }





    @Test
    public void testEqualNumberOfSequences() {

    }



}
