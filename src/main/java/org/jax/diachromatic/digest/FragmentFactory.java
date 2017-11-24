package org.jax.diachromatic.digest;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Class to perform in silico digestion of genome FASTA files. Note that this class borrows some ideas from VPV
 * but the fragments that we create are slightly different compared to VPV because more information is required
 * in the output file.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @version 0.0.1
 */
public class FragmentFactory {
    private static final Logger logger = LogManager.getLogger();
    List<RestrictionEnzyme> restrictionEnzymeList=null;
    List<Fragment> restrictionFragmentList = null;
    List<String> genomeFilePaths = null;

    public String getGenomeDirectoryPath() {
        return genomeDirectoryPath;
    }

    private String genomeDirectoryPath = null;


    public FragmentFactory(String directoryPath) {
        this.genomeDirectoryPath = directoryPath;
        logger.error(String.format("FragmentFactory directory=%s",directoryPath));
        restrictionFragmentList = new ArrayList<>();
        genomeFilePaths = new ArrayList<>();
        identifyFASTAfiles();
        // Note restriction enzyme file is in src/main/resources
        ClassLoader classLoader = FragmentFactory.class.getClassLoader();
        String restrictionEnzymesPath = classLoader.getResource("data/enzymelist.tab").getFile();
        restrictionEnzymeList=RestrictionEnzyme.parseRestrictionEnzymesFromFile(restrictionEnzymesPath);
    }

//    public void createDigest() {
//        for (String path : this.genomeFilePaths) {
//            try {
//                parseFASTAFile(path);
//            } catch (FileNotFoundException e) {
//                e.printStackTrace();
//            }
//        }
//    }

    /** ToDo first iterate over all enzymes and then create restriction fragments like hicup */
//    private void parseFASTAFile(String path) throws FileNotFoundException {
//        IndexedFastaSequenceFile fastaReader = new IndexedFastaSequenceFile(new File(path));
//        String seq = fastaReader.nextSequence().getBaseString();
//
//        Pattern pattern = Pattern.compile(cutpattern);
//        Matcher matcher = pattern.matcher(seq);
//        ArrayList<Integer> cuttingPositionList = new ArrayList<>();
//
//        while (matcher.find()) {
//            int pos = matcher.start() + offset; /* one-based position of first nucleotide after the restriction enzyme cuts */
//            cuttingPositionList.add(pos);
//
//        }
//
//    }


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
            String contigname = null;
            if (fileEntry.isDirectory()) {
                continue;
            } else if (fileEntry.getName().contains("random")) {
                continue; /* skip random contigs! */
            } else if (!fileEntry.getPath().endsWith(".fa")) {
                continue;
            } else {
                this.genomeFilePaths.add(fileEntry.getPath());

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
            logger.fatal(String.format("Could not find enzyme %s, terminating",enzymeName));
            System.exit(1);
        }
    }



    private void cutChromosome(List<RestrictionEnzyme> reList, String chromosomeFilePath) throws Exception {
        IndexedFastaSequenceFile fastaReader = new IndexedFastaSequenceFile(new File(chromosomeFilePath));
        ReferenceSequence refseq = fastaReader.nextSequence();
        Map<String, List<Integer>> cuttingPositionMap=new HashMap<>();
        String seqname = refseq.getName();

        for (RestrictionEnzyme enzyme : reList) {
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
                cuttingPositionList.add(pos);
            }
            cuttingPositionMap.put(enzyme.getPlainSite(), cuttingPositionList); // push array list to map
        }
        // output cuttings for this chromosome to file

    }


}
