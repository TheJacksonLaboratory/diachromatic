package org.jax.diachromatic.io;


import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * The purpose of this class is to encapsulate the FASTA Indexer in a Task.
 * @author Peter Robinson
 * @version 0.2.1 (2017-11-11)
 */
public class FASTAIndexManager {
    private static final Logger logger = LogManager.getLogger();
    /** Path to the directory where we will download and decompress the genome file. */
    List<String> chromosomeFilePaths=null;

    private ProgressBar pbar;


    public FASTAIndexManager(List<String> paths) {
        chromosomeFilePaths=paths;
        pbar = new ProgressBar(0,paths.size());
    }


    public String getContigName(File file) {
        String name=file.getName();
        int i=name.indexOf(".fa");
        if (i>0)
            return name.substring(0,i);
        else
            return name;
    }

    /** Create FAI indices for all FASTA files (only if needed).
     * It is packaged as a Task to allow concurrency*/
    public void indexChromosomes() throws Exception {
        int i=0;
        pbar.print(0);
        for (String path : chromosomeFilePaths) {
            File chromosomeFile = new File(path);
            String contigname=null;
            if (chromosomeFile.isDirectory()) {
                continue;
            } else if (chromosomeFile.getName().contains("random")) {
                continue; /* skip random contigs! */
            } else if (!chromosomeFile.getPath().endsWith(".fa")) {
                continue;
            } else if (fastaFAIalreadyExists(chromosomeFile.getAbsolutePath())) {
                continue;
            } else {
                /* if we get here, we have a FASTA file ending with ".fa" that has not yet been indexed */
                try {
                    FASTAIndexer indexer=new FASTAIndexer(chromosomeFile.getAbsolutePath());
                    pbar.print(++i);
                    indexer.createFASTAindex();
                    indexer.writeFASTAIndex();


                } catch (IOException e) {
                    logger.error("Error encountered while indexing FASTA files");
                    logger.error(e,e);
                    pbar.print(0);
                    throw new Exception(e.getMessage());
                }
                logger.trace("Adding map entry: "+contigname+": "+chromosomeFile.getAbsolutePath());

            }
        }
        pbar.finish();
        return;
    }





    /** Return true if the FASTA index is found already -- no need to repeat! */
    private boolean fastaFAIalreadyExists(String path) {
        File f=new File(String.format("%s.fai",path));
        return f.exists();
    }

}
