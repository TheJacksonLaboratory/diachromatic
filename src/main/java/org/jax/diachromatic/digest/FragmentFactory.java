package org.jax.diachromatic.digest;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
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
 * @version 0.0.1
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


    public FragmentFactory(String directoryPath) {
        this.genomeDirectoryPath = directoryPath;
        logger.error(String.format("FragmentFactory directory=%s",directoryPath));
        restrictionFragmentList = new ArrayList<>();
        genomeFilePaths = new ArrayList<>();
        identifyFASTAfiles();
        // Note restriction enzyme file is in src/main/resources
        restrictionEnzymeList=RestrictionEnzyme.parseRestrictionEnzymes();
    }




    public void digestGenome(List<String> enzymes) throws DiachromaticException {
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
        for (String path : genomeFilePaths) {
            try {
                out = new BufferedWriter(new FileWriter(outfilename));
                out.write(HEADER + "\n");
                cutChromosome(path,out);
                out.close();
            } catch (Exception e) {
                e.printStackTrace();
                logger.fatal("Could not cut chromo. FATAL TODO -- make better exception");
                System.exit(1);
            }
        }
    }


    /**
     * The HTSJDK software requires a FASTA index file for each FASTA file. For instance, for the file {@code example.fa}, HTSJDK
     * would require abn index file that is named {@code example.fa.fai}. We use a
     * {@link FASTAIndexManager} object to go through all the FASTA files in {@link #genomeDirectoryPath} and to
     * generate an index file if needed.
     */
    public void indexFASTAfilesIfNeeded() {
        FASTAIndexManager manager = new FASTAIndexManager(this.genomeFilePaths);
        try {
            System.err.println("indexing fasta files");
            manager.indexChromosomes();
        } catch (Exception e) {
            logger.fatal(String.format("Could not index chromosomes: %s",e.toString()));
            System.exit(1); //TODO make exception
        }
    }

    /** This function will identify FASTA files in the directory {@link #genomeDirectoryPath} by looking for all files
     * with the suffix {@code .fa}. It will add the absolute path of each file on the local file system to the
     * list {@link #genomeFilePaths}.
     */
    private void identifyFASTAfiles() {
        File genomeDir = new File(this.genomeDirectoryPath);
        if (!genomeDir.exists()) {
            logger.error(String.format("Could not find directory \"%s\" with genome FASTA files", this.genomeDirectoryPath));
            System.exit(1); // todo exception
        }
        if (!genomeDir.isDirectory()) {
            logger.error(String.format("%s must be a directory with genome FASTA files", this.genomeDirectoryPath));
            System.exit(1); // todo exception
        }
        for (final File fileEntry : genomeDir.listFiles()) {
            if (fileEntry.isDirectory()) {
                continue;
            } else if (fileEntry.getName().contains("random")) {
                continue; /* skip random contigs! */
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


    /**
     * TODO extend to multiple enzymes
     * TODO throw exception if problems occur
     * @param enzymeName
     */
    public  void cutWithEnzyme(String enzymeName) {
        RestrictionEnzyme re=restrictionEnzymeList.stream().filter(renz -> renz.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (re==null) {
            for (RestrictionEnzyme r: restrictionEnzymeList) {
                System.out.println(r);
            }
            logger.fatal(String.format("Could not find enzyme %s, terminating",enzymeName));
            System.exit(1);
        }
    }



    private void cutChromosome(String chromosomeFilePath, BufferedWriter out) throws Exception {
        IndexedFastaSequenceFile fastaReader = new IndexedFastaSequenceFile(new File(chromosomeFilePath));
        ReferenceSequence refseq = fastaReader.nextSequence();
        ImmutableList.Builder<Fragment> builder = new ImmutableList.Builder<>();
        //List<Fragment> fraglist=new ArrayList<>();
        String seqname = refseq.getName();

        for (Map.Entry<RestrictionEnzyme,Integer> ent : enzyme2number.entrySet()) {
            int enzymeNumber = ent.getValue();
            RestrictionEnzyme enzyme = ent.getKey();
            String cutpat = enzyme.getPlainSite();
            int offset = enzyme.getOffset();

            // note fastaReader refers to one-based numbering scheme.
            String sequence = fastaReader.getSequence(seqname).getBaseString();//(seqname, genomicPos - maxDistToGenomicPosUp, genomicPos + maxDistToGenomicPosDown).getBaseString().toUpperCase();
            Pattern pattern = Pattern.compile(cutpat);
            Matcher matcher = pattern.matcher(sequence);
            ArrayList<Integer> cuttingPositionList = new ArrayList<>();
            /* one-based position of first nucleotide in the entire subsequence returned by fasta reader */
            while (matcher.find()) {
                // replaces matcher.start() - maxDistToGenomicPosUp + offset;
                int pos = matcher.start() + offset; /* one-based position of first nucleotide after the restriction enzyme cuts */
                builder.add(new Fragment(enzymeNumber,pos));
            }
        }
        ImmutableList<Fragment> fraglist = builder.build();
        fraglist = ImmutableList.sortedCopyOf(fraglist);// todo better way of coding this?
        String previousCutEnzyme="None";
        Integer previousCutPosition=1; // start of chromosome
        //Header
        //Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number RE1_Fragment_Number     5'_Restriction_Site     3'_Restriction_Site
        String chromo=seqname;
        int n=0;
        for (Fragment f:fraglist) {
            out.write(String.format("%s\t%d\t%d\t%d\t%s\t%s\n",
                    chromo,
                    previousCutPosition,
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
                previousCutPosition,
                endpos,
                (++n),
                previousCutEnzyme,
                "None"));

    }


}
