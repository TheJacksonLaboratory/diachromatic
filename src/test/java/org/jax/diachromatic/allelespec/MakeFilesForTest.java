package org.jax.diachromatic.allelespec;


import org.junit.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertTrue;

/**
 * This class can be deleted as soon as we are happy with the tests. It is designed to make
 * a "fake" genome and to make fastq files that we can align to make a test BAM file.
 * We run the program like a test but this is just for convenience (ugly hack I guess).
 * Run the test and move the files to some directory. Then...
 * $ bwa index -a bwtsw genome.small.fa
 * Use the following readgroup info
 *
 * @RG\tID:rg1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:H7APS:1:TAAGGCGA $ bwa mem -R'@RG\tID:rg1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:H7APS:1:TAAGGCGA' genome.small.fa sample1.fq sample2.fq >sample.sam
 * $ samtools view -Sb sample.sam > sample.bam
 * $ samtools sort -T sample.sorted sample.bam > sample.sorted.bam
 * $ samtools index sample.sorted.bam
 *
 *
 * This produces a genome sequence like this
 * <pre>
 * GGCCGTCCCT 0
 * TCGCAGATAG 10
 * TGGTACTCTG 20
 * CCTTGCCCCC 30
 * CTATCCAATC 40
 * AAGCTTCGGC 50
 * GCTTGCGTAT 60
 * GCCGGACAAC 70
 * GCTCGTTGTA 80
 * TCTCCCCCAG 90
 * ATAAAGCTTC 100
 * ATGTGGAAGC 110
 * TTGAGTGCTT 120
 * CTCGCTGGCC 130
 * GTCTGGTCTG 140
 * ATGTTTTAGG 150
 * AGAAGCTTGC 160
 * GCAAAATCTC 170
 * CAGACGGGGT 190
 * TCCGCTGA   200
 * </pre>
 * There are 19 reads that overlap with position 20. The wildtype sequence
 * at position 20 is a 'T'. 3 of the reads have an 'A' at this position. There
 * are 19 reads that overlap with position 70. The wildtype sequence at this
 * position is a 'G'; 7 of the reads have a 'C' at this position. This can be
 * visualized with
 * $ samtools tview sample.sorted.bam genome.small.fa
 *
 */
public class MakeFilesForTest {

    private final static int READLEN = 20;

    private List<String> fastq1;
    private List<String> fastq2;


    private char makeRandomNT(double r) {
        if (r < 0.25) {
            return 'A';
        } else if (r < 0.5) {
            return 'C';
        } else if (r < 0.75) {
            return 'G';
        } else {
            return 'T';
        }
    }

    private String makeran(Random random, int len) {
        char[] chars = new char[len];
        for (int i = 0; i < len; i++) {
            char c = makeRandomNT(random.nextDouble());
            chars[i] = c;
        }
        String seq = new String(chars);
        return seq;
    }


    @Test
    public void createGenomeFile() {
        long seed = 42L;
        Random rand = new Random(seed);
        String hindIII = "AAGCTT";
        StringBuilder sb = new StringBuilder();
        sb.append(makeran(rand, 50));
        sb.append(hindIII);
        sb.append(makeran(rand, 60));
        sb.append(hindIII);
        sb.append(makeran(rand, 40));
        sb.append(hindIII);
        sb.append(makeran(rand, 30));
        String genome = sb.toString();
        System.out.println(genome);
        assertTrue(genome.length() > 10);
        int found = 0;
        Pattern pat = Pattern.compile(hindIII);
        Matcher matcher = pat.matcher(genome);
        List<Integer> poslist = new ArrayList<>();
        while (matcher.find()) {
            found++;
            int pos = matcher.start();
            poslist.add(pos);
        }
        assertEquals(4, found);
        // output genome file
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("genome.small.fa"));
            writer.write(">chrZ\n");
            writer.write(genome);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        // output some paired reads
        // we will go for fragments 1-3 and fragments 2-4
        fastq1 = new ArrayList<>();
        fastq2 = new ArrayList<>();
        output1and3(genome, poslist);
        output2and4(genome, poslist);
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("sample1.fq"));
            for (String line : fastq1) {
                writer.write(line + "\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter("sample2.fq"));
            for (String line : fastq2) {
                writer.write(line + "\n");
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("[INFO] Done writing genome and two fastq files");
    }


    private char getOtherNT(char c) {
        switch (c) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
            default:
                return 'N';// should never happen
        }
    }

    private char getQual(Random rand) {
        final String qualstring = "!\"#$%&'()*+,-./0123456789:;<=>?ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz";
        final int len = qualstring.length();
        int i = rand.nextInt(len);
        return qualstring.charAt(i);
    }


    private String getReadName(boolean forward, int x, int y) {
        String[] fields = {"HWI-D00119", "50", "H7AP8ADXY", "1", "42", "2100", "2300", "2", "N", "0", "TAAGGCGA"};
        //if (forward) fields[7]="1";
        fields[5] = String.valueOf(x);
        fields[6] = String.valueOf(y);
        return String.format("@%s", String.join(":", fields));
    }

    private String getQualString(int len, Random rand) {
        char[] qual = new char[len];
        for (int i = 0; i < len; i++) {
            qual[i] = getQual(rand);
        }
        return new String(qual);
    }


    private void output1and3(String genome, List<Integer> hindIIIpos) {
        // Add a mutation at position 20 and have reads that start at position 12-20
        // three each
        // chance of a mutation is 0.2
        char mutchar = getOtherNT(genome.charAt(20));
        String mutant = genome.substring(0, 20) + mutchar + genome.substring(21);
        System.err.println(genome);
        System.err.println(mutant);
        double MUTATION_THRESHOLD = 0.2;
        Random rand = new Random(13);
        int x = 100;
        int y = 50;
        for (int i = 2; i < 21; i++) {
            for (int j = 0; j < 1; j++) {
                int endpos = i + READLEN;
                double r = rand.nextDouble();
                String readname = getReadName(true, x, y);
                String qual = getQualString(READLEN, rand);
                String subseq;
                if (r < MUTATION_THRESHOLD) {
                    subseq = mutant.substring(i, endpos);
                } else {
                    subseq = genome.substring(i, endpos);
                }
                this.fastq1.add(readname);
                this.fastq1.add(subseq);
                this.fastq1.add("+");
                this.fastq1.add(qual);
                x++;
            }
        }
        // the reads in the thrid fragment do not have a mutation
        x = 100; // reset so the numbers are the same
        for (int i = 110; i < 129; i++) {
            String readname = getReadName(false, x, y);
            String qual = getQualString(READLEN, rand);
            int endpos = i + READLEN;
            String subseq = genome.substring(i, endpos);
            this.fastq2.add(readname);
            this.fastq2.add(subseq);
            this.fastq2.add("+");
            this.fastq2.add(qual);
            x++;
        }
    }


    private void output2and4(String genome, List<Integer> hindIIIpos) {
        // Add a mutation at position 20 and have reads that start at position 12-20
        // three each
        // chance of a mutation is 0.5
        char mutchar = getOtherNT(genome.charAt(70));
        String mutant = genome.substring(0, 70) + mutchar + genome.substring(71);
        System.err.println(genome);
        System.err.println(mutant);
        double MUTATION_THRESHOLD = 0.5;
        Random rand = new Random(13);
        int x = 100;
        int y = 50;
        for (int i = 52; i < 71; i++) {
            int endpos = i + READLEN;
            double r = rand.nextDouble();
            String readname = getReadName(true, x, y);
            String qual = getQualString(READLEN, rand);
            String subseq;
            if (r < MUTATION_THRESHOLD) {
                System.out.println("MUTANT");
                subseq = mutant.substring(i, endpos);
            } else {
                subseq = genome.substring(i, endpos);
            }
            this.fastq1.add(readname);
            this.fastq1.add(subseq);
            this.fastq1.add("+");
            this.fastq1.add(qual);
            x++;

        }
        // the reads in the fourth fragment do not have a mutation
        x = 100;
        for (int i = 155; i < 174; i++) {
            for (int j = 0; j < 1; j++) {
                String readname = getReadName(true, x, y);
                String qual = getQualString(READLEN, rand);
                int endpos = i + READLEN;
                String subseq = genome.substring(i, endpos);
                this.fastq2.add(readname);
                this.fastq2.add(subseq);
                this.fastq2.add("+");
                this.fastq2.add(qual);
                x++;
            }
        }


    }


}
