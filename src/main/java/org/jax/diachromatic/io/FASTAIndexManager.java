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
 * @version 0.2.2 (2018-03-19)
 */
public class FASTAIndexManager {
    private static final Logger logger = LogManager.getLogger();
    /** Path to the directory where we will download and decompress the genome file. */
    private List<String> chromosomeFilePaths=null;


    public FASTAIndexManager(List<String> paths) {
        chromosomeFilePaths=paths;
    }

    public static void indexChromosome(String path) throws Exception {
        if (fastaFAIalreadyExists(path)) {
            logger.trace("Cowardly refusing to index "+path+" because FAI file was found");
            return;
        }
        try {
            logger.trace("Indexing chromosome "+path);
            FASTAIndexer indexer=new FASTAIndexer(path);
            indexer.createFASTAindex();
            indexer.writeFASTAIndex();
        } catch (IOException e) {
            logger.error("Error encountered while indexing FASTA files");
            logger.error(e,e);
            throw new Exception(e.getMessage());
        }
    }

    /** Return true if the FASTA index is found already -- no need to repeat! */
    private static boolean fastaFAIalreadyExists(String path) {
        File f=new File(String.format("%s.fai",path));
        return f.exists();
    }

}
