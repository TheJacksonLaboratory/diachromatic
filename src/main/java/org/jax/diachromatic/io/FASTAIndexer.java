package org.jax.diachromatic.io;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.RandomAccessFile;

/**
 * This class creates a FASTA index for the FASTA file passed to the constructor.
 * Note that in the future it would be better to use HTSJDK, but the FASTA Indexing function
 * is not available with the current version in maven central.
 * * The fasta.fai is the fasta index,.
 * <ol><li>Column 1: The contig name.</li>
 * <li>Column 2: The number of bases in the contig</li>
 * <li>Column 3: The byte index of the file where the contig sequence begins.</li>
 * <li>Column 4: bases per line in the FASTA file</li>
 * <li>Column 5: bytes per line in the FASTA file</li>
 * </ol>
 * For instance, the fai for human chr15 is
 *  <pre>
 *     chr15	102531392	7	50	51
 *  </pre>
 * Note that for the UCSC files, there is one sequence per file, and thus the fai files have only one line.
 * Also note that I did not use the HTSJDK indexer because it is not available in the versions of the library in maven
 * central. We can later use their implementation if this changes. The class should create a fasta index equivalent to
 * what would be produced by <pre>$ samtools faidx test.fa</pre>.
 * @author Peter Robinson
 * @version 0.0.2 (2017-07-24).
 */
public class FASTAIndexer  {
    private static final Logger logger = LogManager.getLogger();
    /** Path to the directory where we will download and decompress the genome file. */
    private String fastaPath=null;
    private String fastaFaiPath=null;
    private String contigname;

    private long n_bases;
    private long byte_index;
    private long bases_per_line;
    private long bytes_per_line;

    /**
     *
     * @param path Path to the FASTA file to be indexed.
     */
    public FASTAIndexer(String path) {
        this.fastaPath=path;
        this.fastaFaiPath=String.format("%s.fai",path);
        logger.trace("Indexing fasta file at directory "+path);
    }


    public String getContigname(){return contigname;}
    public long getN_bases() {
        return n_bases;
    }

    public long getByte_index() {
        return byte_index;
    }

    public long getBases_per_line() {
        return bases_per_line;
    }

    public long getBytes_per_line() {
        return bytes_per_line;
    }


    /**
    * Extract data for a FASTA index from the file passed to the function
     * We return the contig name under the assumption that the contig name
     * matches the file name exactly except for the .fa suffix.
     */
    public void createFASTAindex() throws IOException {
        RandomAccessFile file = new RandomAccessFile(this.fastaPath, "r");
        file.seek(0);
        String header = file.readLine();
        if (!header.startsWith(">")) {
            logger.error("FASTA header line did not start with >: " + header);
            return;
        }
        this.contigname = header.substring(1);
        this.byte_index = file.getFilePointer(); /* this is the offset in bytes of the first sequence line. */
        String line = file.readLine();
        this.bases_per_line = line.length();
        long offset = file.getFilePointer();
        this.bytes_per_line = offset - this.byte_index;
        this.n_bases = bases_per_line;
        while ((line = file.readLine()) != null) {
            this.n_bases += line.length();
        }
        file.close();
    }

    public void writeFASTAIndex()  throws IOException{
        BufferedWriter bw = new BufferedWriter(new FileWriter(this.fastaFaiPath));
        bw.write(String.format("%s\t%d\t%d\t%d\t%d\n",contigname,n_bases,byte_index,bases_per_line,bytes_per_line));
        bw.close();
    }


}
