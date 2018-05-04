package org.jax.diachromatic.command;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.Bowtie2Runner;
import org.jax.diachromatic.align.Digest;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
/**
 * Class to coordinate bowtie2-mapping of the truncated FASTQ files followed by Q/C and filtering of the
 * mapped reads.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */
public class AlignCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to the bowite2 exectuable, e.g., {@code /usr/bin/bowtie2}. */
    private final String bowtiepath;

    /** Path to the bowtie2 index files. Note that the index is made up of multiple files, e.g.,
     * hg19.1.bt2,  hg19.3.bt2,  hg19.rev.1.bt2, hg19.2.bt2,  hg19.4.bt2,  hg19.rev.2.bt2. Assuming all files
     * are in a directory called {@code /path/to/index/}, this parameter should be {@code /path/to/index/hg19}.*/
    private final String pathToBowtieIndex;

    /** Path to the forward truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    private String pathToInputFastq1 = null;

    /** Path to the reverse truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    private String pathToInputFastq2 = null;


    /** Path to the genome digest file produced by {@link org.jax.diachromatic.command.DigestCommand}.*/
    private String digestFile = null;

    /** Path to BED file containing the coordinates of active digests {@link org.jax.diachromatic.command.DigestCommand}.*/
    private String activeDigestsFile = null;

    //** if this is set, an extra BAM file containg the rejected read pairs will be created */
    private final boolean outputRejectedReads;

    private String outputPathPrefix = null;

    public AlignCommand(String bowtie, String btIndexPath, String inputFastqPath1, String inputFastqPath2, String digest, String activeDigests,
                        boolean outputRejected, String outputPathPrefix) {
        this.bowtiepath =bowtie;
        pathToBowtieIndex=btIndexPath;
        pathToInputFastq1 =inputFastqPath1;
        pathToInputFastq2 =inputFastqPath2;
        digestFile=digest;
        activeDigestsFile=activeDigests;
        outputRejectedReads=outputRejected;
        this.outputPathPrefix=outputPathPrefix;
    }

    public void execute() {
        String samFile1="tempoutname1.sam";
        String samFile2="tempoutname2.sam";
        logger.trace(String.format("About to read digests from %s",digestFile));
        Map<String,List<Digest>> digestmap = Digest.readDigests(digestFile, activeDigestsFile);
        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq1,samFile1);
            runner.run();
            Bowtie2Runner runner2 = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq2,samFile2);
            runner2.run();
            Aligner pairer = new Aligner(samFile1,samFile2,digestmap,outputRejectedReads,outputPathPrefix);
            pairer.inputSAMfiles();
            pairer.printStatistics();
            File file = new File(samFile1);
            file.delete();
            file = new File(samFile2);
            file.delete();
        } catch (DiachromaticException | IOException e){
            e.printStackTrace();
        }

    }
    @Override
    public String toString() {return "diachromatic:align";} //???
}