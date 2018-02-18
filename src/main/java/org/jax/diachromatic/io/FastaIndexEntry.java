package org.jax.diachromatic.io;


import org.jax.diachromatic.exception.UnindexableFastaFileException;

/** This class represents an entry of the fasta index (i.e. a line of the .fai file). It has been adapted
 * from classes from the ASCII Genome The MIT License (MIT); Copyright (c) 2016 Dario Beraldi
 *
 * This class creates a FASTA index for the FASTA file passed to the constructor.
 * Note that in the future it would be better to use HTSJDK, but the FASTA Indexing function
 * is not available with the current version in maven central.
 * * The fasta.fai is the fasta index,.
 * <ol>
 * <li>Column 1: The contig name.</li>
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
public class FastaIndexEntry {
    private String seqName= null;
    protected int seqLength= 0;
    protected long byteOffset= 0; // Byte position where the sequence starts.
    protected int lineLength= 0; // Only nucleotides
    protected int lineFullLength= 0; // Including line terminators.

    protected void makeSeqNameFromRawLine(String line) throws UnindexableFastaFileException {
        if( ! line.startsWith(">") ){
            throw new UnindexableFastaFileException("Invalid name: Does not start with '>'");
        }
        this.seqName = line.substring(1).trim().replaceAll("\\s.*", "");
    }

    protected String getSeqName(){
        return seqName;
    }

    @Override
    public String toString(){
        return seqName + "\t" + seqLength + "\t" + byteOffset + "\t" + lineLength + "\t" + lineFullLength;
    }
}
