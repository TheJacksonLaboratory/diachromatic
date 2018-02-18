package org.jax.diachromatic.digest;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.Diachromatic;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.io.FASTAIndexManager;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to perform in silico digestion of genome FASTA files. Note that this class borrows some ideas from VPV
 * but the fragments that we create are slightly different compared to VPV because more information is required
 * in the output file.
 * <pre>
 * Genome:testgenome       Restriction_Enzyme1:BgIII [A^GATCT]     Restriction_Enzyme2:None        Hicup digester version 0.5.10
 * Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number RE1_Fragment_Number     5'_Restriction_Site     3'_Restriction_Site
 * chrUn_KI270745v1        1       1861    1       1       None    Re1
 * chrUn_KI270745v1        1862    29661   2       2       Re1     Re1
 * chrUn_KI270745v1        29662   35435   3       3       Re1     Re1
 * chrUn_KI270745v1        35436   40296   4       4       Re1     Re1
 * chrUn_KI270745v1        40297   41891   5       5       Re1     None
 * </pre>
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @version 0.0.3
 */
public class FragmentFactory {
    private static final Logger logger = LogManager.getLogger();
    private List<RestrictionEnzyme> restrictionEnzymeList = null;
    private Map<Integer, RestrictionEnzyme> number2enzyme = null;
    private Map<RestrictionEnzyme, Integer> enzyme2number = null;
    private List<Fragment> restrictionFragmentList = null;
    /**
     * File handle for the output of the restriction fragments.
     */
    private BufferedWriter out = null;
    /**
     * Name of output file. TODO set this dynamically
     */
    private String outfilename = "hicupCloneDigest.txt";
    /**
     * A counter used to output each 1000 fragment as we go as feedback to the user.
     */
    private int counter = 1;

    public String getGenomeDirectoryPath() {
        return genomeFastaFilePath;
    }

    /**
     * Path to the genome FASTA file, e.g., hg19.fa. Note that expect there to be a corresponding index, e.g., hg19.fa.fai,
     * or we will build one.
     */
    private final String genomeFastaFilePath;
    /**
     * Header of the output file.
     */
    private static final String HEADER = "Chromosome\tFragment_Start_Position\t" +
            "Fragment_End_Position\tFragment_Number\t5'_Restriction_Site\t3'_Restriction_Site";

    public FragmentFactory(String fastFilePath, String outfile) {
        this.genomeFastaFilePath = fastFilePath;
        outfilename = outfile;
        logger.trace("FragmentFactory will create digest from fasta file={}", fastFilePath);
        restrictionFragmentList = new ArrayList<>();
        // Note restriction enzyme file is in src/main/resources
        restrictionEnzymeList = RestrictionEnzyme.parseRestrictionEnzymes();
    }


    /**
     * Produces a digest file starting from a single genome fasta file, e.g., hg19.fa.
     *
     * @param enzymes
     * @throws DiachromaticException
     */
    public void digestGenome(List<String> enzymes) throws DiachromaticException {
        number2enzyme = new HashMap<>();
        enzyme2number = new HashMap<>();
        int n = 0;
        for (String enzym : enzymes) {
            RestrictionEnzyme re = restrictionEnzymeList.stream().
                    filter(x -> enzym.equalsIgnoreCase(x.getName())).
                    findFirst().orElse(null);
            if (re == null) {
                throw new DiachromaticException(String.format("Did not recognize restriction enzyme \"%s\"", enzym));
            } else {
                n++;
                number2enzyme.put(n, re);
                enzyme2number.put(re, n);
            }
        }
        try {
            out = new BufferedWriter(new FileWriter(outfilename));
            out.write(HEADER + "\n");
            FASTAIndexManager.indexChromosome(this.genomeFastaFilePath);
            cutGenome(this.genomeFastaFilePath, out);
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new DiachromaticException(String.format("Could not digest chromosomes: %s", e.toString()));
        }

    }


    /**
     * Cut an entire genome file.
     *
     * @param genomeFasta Path to the genome fasta file, e.g., hg19.fa
     * @param out         file handle to which we will write the digest file
     * @throws DiachromaticException if we cannot write the digest file
     */

    private void cutGenome(String genomeFasta, BufferedWriter out) throws DiachromaticException {
        logger.trace(String.format("cutting chromosomes %s", genomeFasta));
        IndexedFastaSequenceFile fastaReader = null;
        try {
            fastaReader = new IndexedFastaSequenceFile(new File(genomeFasta));
        } catch (Exception e) {
            throw new DiachromaticException(String.format("Could not create FAI file for %s [%s]", genomeFasta, e.toString()));
        }
        ReferenceSequence refseq;
        try {
            while ((refseq = fastaReader.nextSequence()) != null) {
                String seqname = refseq.getName();
                String sequence = fastaReader.getSequence(seqname).getBaseString().toUpperCase();
                cutChromosome(seqname, sequence, out);
            }
        } catch (IOException e) {
            throw new DiachromaticException("Unable to write the digest file: " + e.getMessage());
        }
    }

    /**
     * Write the corresponding part of the digest file for one chromosome
     * @param seqname Name of the chromosome
     * @param sequence Sequence of the chromosome
     * @param out file handle to which we will write the digest file
     * @throws DiachromaticException
     * @throws IOException
     */
    private void cutChromosome(String seqname, String sequence, BufferedWriter out) throws DiachromaticException,IOException {
        ImmutableList.Builder<Fragment> builder = new ImmutableList.Builder<>();
        for (Map.Entry<RestrictionEnzyme, Integer> ent : enzyme2number.entrySet()) {
            int enzymeNumber = ent.getValue();
            RestrictionEnzyme enzyme = ent.getKey();
            String cutpat = enzyme.getPlainSite();
            int offset = enzyme.getOffset();
            // note fastaReader refers to one-based numbering scheme.
            Pattern pattern = Pattern.compile(cutpat);
            Matcher matcher = pattern.matcher(sequence);
            /* one-based position of first nucleotide in the entire subsequence returned by fasta reader */
            while (matcher.find()) {
                int pos = matcher.start() + offset; /* one-based position of first nucleotide after the restriction enzyme cuts */
                if (++counter % 1000 == 0) {
                    logger.trace(String.format("Added %d-th fragment", counter));
                }
                builder.add(new Fragment(enzymeNumber, pos));
            }
        }
        ImmutableList<Fragment> fraglist = builder.build();
        fraglist = ImmutableList.sortedCopyOf(fraglist);// todo better way of coding this?
        String previousCutEnzyme = "None";
        Integer previousCutPosition = 0; // start of chromosome
        //Header
        //Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number RE1_Fragment_Number     5'_Restriction_Site     3'_Restriction_Site
        String chromo = seqname;
        int n = 0;
        for (Fragment f : fraglist) {
            out.write(String.format("%s\t%d\t%d\t%d\t%s\t%s\n",
                    chromo,
                    (previousCutPosition + 1),
                    f.position,
                    (++n),
                    previousCutEnzyme,
                    number2enzyme.get(f.enzymeNumber).getName()));
            previousCutEnzyme = number2enzyme.get(f.enzymeNumber).getName();
            previousCutPosition = f.position;
        }
        // output last fragment also
        // No cut ("None") at end of chromosome
        int endpos = seqname.length();
        out.write(String.format("%s\t%d\t%d\t%d\t%s\t%s\n",
                chromo,
                (previousCutPosition + 1),
                endpos,
                (++n),
                previousCutEnzyme,
                "None"));

    }


}
