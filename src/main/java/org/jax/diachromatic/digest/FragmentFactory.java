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
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to perform in silico digestion of genome FASTA files. Note that this class borrows some ideas from VPV
 * but the fragments that we create are slightly different compared to VPV because more information is required
 * in the output file.
 * TODO decide on format. The following is the format of HiCUP. TODO question is whether we want to support the
 * Re2 (used instead of sonication)? Maybe not, because we are going to be using 4-cutters with smaller fragments anyway
 * and thus it is pretty unlikely that anybody will be using a second enzyme to fragment our fragments.
 *
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
 * @version 0.0.2
 */
public class FragmentFactory {
    private static final Logger logger = LogManager.getLogger();
    List<RestrictionEnzyme> restrictionEnzymeList=null;
    Map<Integer,RestrictionEnzyme> number2enzyme =null;
    Map<RestrictionEnzyme,Integer> enzyme2number=null;
    List<Fragment> restrictionFragmentList = null;
    List<String> genomeFilePaths = null;
    /** File handle for the output of the restriction fragments. */
    private BufferedWriter out = null;
    /** Name of output file. TODO set this dynamically */
    private String outfilename="hicupCloneDigest.txt";

    public String getGenomeDirectoryPath() {
        return genomeDirectoryPath;
    }

    private String genomeDirectoryPath = null;
    /** Header of the output file. */
    private static final String HEADER="Chromosome\tFragment_Start_Position\t" +
            "Fragment_End_Position\tFragment_Number\t5'_Restriction_Site\t3'_Restriction_Site";


    public FragmentFactory(String directoryPath, String outfile) {
        this.genomeDirectoryPath = directoryPath;
        outfilename=outfile;
        logger.trace(String.format("FragmentFactory directory=%s",directoryPath));
        restrictionFragmentList = new ArrayList<>();
        genomeFilePaths = new ArrayList<>();
        // Note restriction enzyme file is in src/main/resources
        restrictionEnzymeList=RestrictionEnzyme.parseRestrictionEnzymes();
    }




    public void digestGenome(List<String> enzymes) throws DiachromaticException {
        identifyFASTAfiles();
        number2enzyme =new HashMap<>();
        enzyme2number=new HashMap<>();
        int n=0;
        for (String enzym : enzymes) {
            RestrictionEnzyme re = restrictionEnzymeList.stream().
                    filter( x ->  enzym.equalsIgnoreCase(x.getName()) ).
                    findFirst().orElse(null);
            if (re==null) {
                throw new DiachromaticException(String.format("Did not recognize restriction enzyme \"%s\"",enzym));
            } else {
                n++;
                number2enzyme.put(n,re);
                enzyme2number.put(re,n);
            }
        }
        try {
            out = new BufferedWriter(new FileWriter(outfilename));
            out.write(HEADER + "\n");
            for (String path : genomeFilePaths) {
                FASTAIndexManager.indexChromosome(path);
                cutChromosome(path, out);
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
            throw new DiachromaticException(String.format("Could not digest chromosomes: %s", e.toString()));
        }

    }



    /** This function will identify FASTA files in the directory {@link #genomeDirectoryPath} by looking for all files
     * with the suffix {@code .fa}. It will add the absolute path of each file on the local file system to the
     * list {@link #genomeFilePaths}.
     */
    private void identifyFASTAfiles() throws DiachromaticException {
        File genomeDir = new File(this.genomeDirectoryPath);
        if (!genomeDir.exists()) {
            throw new DiachromaticException(String.format("Could not find directory \"%s\" with genome FASTA files", this.genomeDirectoryPath));
        }
        if (!genomeDir.isDirectory()) {
            throw new DiachromaticException(String.format("%s must be a directory with genome FASTA files", this.genomeDirectoryPath));
        }
        for (final File fileEntry : genomeDir.listFiles()) {
            if (fileEntry.isDirectory()) {
                continue;
            } else if (!fileEntry.getPath().endsWith(".fa")) {
                continue;
            } else {
                this.genomeFilePaths.add(fileEntry.getAbsolutePath());
            }
        }
    }


    public int getGenomeFileCount() {
        return genomeFilePaths.size();
    }

    int counter=1;

    private void cutChromosome(String chromosomeFilePath, BufferedWriter out) throws Exception {
        logger.trace(String.format("cutting chromosomes %s",chromosomeFilePath ));
        IndexedFastaSequenceFile fastaReader=null;
        try {


             fastaReader = new IndexedFastaSequenceFile(new File(chromosomeFilePath));
        } catch (Exception e) {
            throw  new DiachromaticException(String.format("Could not create FAI file for %s [%s]",chromosomeFilePath,e.toString()));
        }
        ReferenceSequence refseq = fastaReader.nextSequence();
        ImmutableList.Builder<Fragment> builder = new ImmutableList.Builder<>();
        String seqname = refseq.getName();

        for (Map.Entry<RestrictionEnzyme,Integer> ent : enzyme2number.entrySet()) {
            int enzymeNumber = ent.getValue();
            RestrictionEnzyme enzyme = ent.getKey();
            String cutpat = enzyme.getPlainSite();
            int offset = enzyme.getOffset();

            // note fastaReader refers to one-based numbering scheme.
            String sequence = fastaReader.getSequence(seqname).getBaseString().toUpperCase();//(seqname, genomicPos - maxDistToGenomicPosUp, genomicPos + maxDistToGenomicPosDown).getBaseString().toUpperCase();
            Pattern pattern = Pattern.compile(cutpat);
            Matcher matcher = pattern.matcher(sequence);
           // ArrayList<Integer> cuttingPositionList = new ArrayList<>();
            /* one-based position of first nucleotide in the entire subsequence returned by fasta reader */
            while (matcher.find()) {
                // replaces matcher.start() - maxDistToGenomicPosUp + offset;
                int pos = matcher.start() + offset; /* one-based position of first nucleotide after the restriction enzyme cuts */
               // logger.trace(String.format("Adding %d to search for %s",pos,cutpat));
                if (counter%1000==0) {
                    System.out.println(String.format("Added %d th fragment",counter ));
                }
                builder.add(new Fragment(enzymeNumber,pos));
            }
        }
        ImmutableList<Fragment> fraglist = builder.build();
        fraglist = ImmutableList.sortedCopyOf(fraglist);// todo better way of coding this?
        String previousCutEnzyme="None";
        Integer previousCutPosition=0; // start of chromosome
        //Header
        //Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number RE1_Fragment_Number     5'_Restriction_Site     3'_Restriction_Site
        String chromo=seqname;
        int n=0;
        for (Fragment f:fraglist) {
            out.write(String.format("%s\t%d\t%d\t%d\t%s\t%s\n",
                    chromo,
                    (previousCutPosition+1),
                    f.position,
                    (++n),
                    previousCutEnzyme,
                    number2enzyme.get(f.enzymeNumber).getName()));
            previousCutEnzyme=number2enzyme.get(f.enzymeNumber).getName();
            previousCutPosition=f.position;
        }
        // output last fragment also
        // No cut ("None") at end of chromosome
        int endpos = refseq.length();
        out.write(String.format("%s\t%d\t%d\t%d\t%s\t%s\n",
                chromo,
                (previousCutPosition+1),
                endpos,
                (++n),
                previousCutEnzyme,
                "None"));

    }


}
