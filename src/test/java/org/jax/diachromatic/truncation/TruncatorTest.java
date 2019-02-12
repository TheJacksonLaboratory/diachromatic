package org.jax.diachromatic.truncation;

import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.util.Pair;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;


import static org.junit.jupiter.api.Assertions.assertEquals;

<<<<<<< HEAD
import static org.junit.Assert.assertEquals;
=======
>>>>>>> master

public class TruncatorTest {

    private static FastqPairParser parser=null;
    private static String fastq_1;
    private static String fastq_2;
    private static String ligationSequence;

    @BeforeAll
    public  static void init() {
      ClassLoader classLoader = TruncatorTest.class.getClassLoader();
      fastq_1 = classLoader.getResource("data/fastq/test1.fastq").getFile();
      fastq_2 = classLoader.getResource("data/fastq/test2.fastq").getFile();
      RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
      ligationSequence=Truncator.fillEnd(hindIII);
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
    public void testBmtI () {
        //BmtI	GCTAG^C
        // G-CTAG-CTAG-C = GCTAGCTAGC
        RestrictionEnzyme bmtI = new RestrictionEnzyme("BmtI","GCTAG^C" );
        String ligationSequence = Truncator.fillEnd(bmtI);
        String expected = "GCTAGCTAGC";
        assertEquals(expected,ligationSequence);
    }


    @Test
    public void testFaeI () {
        //FaeI	CATG^
        // CATG-CATG = CATGCATG
        RestrictionEnzyme faeI = new RestrictionEnzyme("FaeI","CATG^" );
        String ligationSequence = Truncator.fillEnd(faeI);
        String expected = "CATGCATG";
        assertEquals(expected,ligationSequence);
    }


    @Test
    public void testPvuI () {
        // PvuI	CGAT^CG
        // CG-AT-AT-CG = CGATATCG
        RestrictionEnzyme pvuI = new RestrictionEnzyme("PvuI","CGAT^CG" );
        String ligationSequence = Truncator.fillEnd(pvuI);
        String expected = "CGATATCG";
        assertEquals(expected,ligationSequence);
    }

    /*
    Search in test1.fastq reveals three instances of AAGCTAGCTT. Thus,
    read1 should be truncated 3 times. There are no hits for this sequence
    int test2.fastq
     */
    @Test
    public void testTruncationCount() throws DiachromaticException{
        RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
        String ligationSequence = Truncator.fillEnd(hindIII);
        PotentiallyTruncatedFastQRecord.setLigationSequence(ligationSequence);
        PotentiallyTruncatedFastQRecord.setRestrictionSequence(hindIII.getPlainSite());
        PotentiallyTruncatedFastQRecord.setDanglingSequence(hindIII.getDanglingEndSequence());
        parser=new FastqPairParser(fastq_1,fastq_2,ligationSequence);
        while (parser.hasNextPair()) {
            Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> pair = parser.getNextPair();
            PotentiallyTruncatedFastQRecord fqr1=pair.first;
            PotentiallyTruncatedFastQRecord fqr2=pair.second;
        }
        assertEquals(3,parser.getReadOneTruncated());
        assertEquals(0,parser.getReadTwoTruncated());
    }


    /*
       Search in test1.fastq reveals three instances of AAGCTAGCTT. Thus,
       read1 should be truncated 3 times. There are no hits for this sequence
       int test2.fastq
        */
    @Test
    public void testReadsProcessed() throws DiachromaticException {
        RestrictionEnzyme hindIII = new RestrictionEnzyme("HindIII","A^AGCTT");
        String ligationSequence = Truncator.fillEnd(hindIII);
        PotentiallyTruncatedFastQRecord.setLigationSequence(ligationSequence);
        PotentiallyTruncatedFastQRecord.setRestrictionSequence(hindIII.getPlainSite());
        PotentiallyTruncatedFastQRecord.setDanglingSequence(hindIII.getDanglingEndSequence());
        parser=new FastqPairParser(fastq_1,fastq_2,ligationSequence);
        while (parser.hasNextPair()) {
            Pair<PotentiallyTruncatedFastQRecord, PotentiallyTruncatedFastQRecord> pair = parser.getNextPair();
            PotentiallyTruncatedFastQRecord fqr1=pair.first;
            PotentiallyTruncatedFastQRecord fqr2=pair.second;
        }
        assertEquals(12,parser.getnReadsProcessed());
    }

}
