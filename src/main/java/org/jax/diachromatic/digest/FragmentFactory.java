package org.jax.diachromatic.digest;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.apache.log4j.Logger;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
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
    static Logger logger = Logger.getLogger(FragmentFactory.class.getName());
    List<Fragment> restrictionFragmentList = null;
    List<String> genomeFilePaths = null;

    public String getGenomeDirectoryPath() {
        return genomeDirectoryPath;
    }

    private String genomeDirectoryPath = null;
    /**
     * Todo make setable and make an array
     */
    private String cutpattern = "AGTC";
    /**
     * Todo make setable and make an array
     */
    private int offset = 0;


    public FragmentFactory(String directoryPath) {
        this.genomeDirectoryPath = directoryPath;
        restrictionFragmentList = new ArrayList<>();
        genomeFilePaths = new ArrayList<>();
        identifyFASTAfiles();
    }

    public void createDigest() {
        for (String path : this.genomeFilePaths) {
            try {
                parseFASTAFile(path);
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }
    }

    /** ToDo first iterate over all enzymes and then create restriction fragments like hicup */
    private void parseFASTAFile(String path) throws FileNotFoundException {
        IndexedFastaSequenceFile fastaReader = new IndexedFastaSequenceFile(new File(path));
        String seq = fastaReader.nextSequence().getBaseString();

        Pattern pattern = Pattern.compile(cutpattern);
        Matcher matcher = pattern.matcher(seq);
        ArrayList<Integer> cuttingPositionList = new ArrayList<>();

        while (matcher.find()) {
            int pos = matcher.start() + offset; /* one-based position of first nucleotide after the restriction enzyme cuts */
            cuttingPositionList.add(pos);

        }

    }


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


}
