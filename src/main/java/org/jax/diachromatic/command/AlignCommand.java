package org.jax.diachromatic.command;

import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.Bowtie2Runner;

import java.io.File;
import java.io.IOException;
import java.util.Random;

/**
 * Class to coordinate bowtie2-mapping of the truncated FASTQ files followed by Q/C and filtering of the
 * mapped reads.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */
@Parameters(commandDescription = "Align Hi-C or capture Hi-C reads with bowtie2")
public class AlignCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to the bowtie2 executable, e.g., {@code /usr/bin/bowtie2}. */
    @Parameter(names={"-b","--bowtie-path"},required = true, description ="path to bowtie2" )
    private String bowtiepath;

    /** Path to the bowtie2 index files. Note that the index is made up of multiple files, e.g.,
     * hg19.1.bt2,  hg19.3.bt2,  hg19.rev.1.bt2, hg19.2.bt2,  hg19.4.bt2,  hg19.rev.2.bt2. Assuming all files
     * are in a directory called {@code /path/to/index/}, this parameter should be {@code /path/to/index/hg19}.*/

    @Parameter(names={"-i", "--bowtie-index"}, required = true, description ="path to bowtie2 index")
    private String pathToBowtieIndex;

    /** Path to the forward truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    @Parameter(names={"-q","fastq-r1"}, required = true,description = "path to forward FASTQ input file")
    private String pathToInputFastq1 = null;

    /** Path to the reverse truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    @Parameter(names={"-r","fastq-r2"}, required = true,description = "path to reverse FASTQ input file")
    private String pathToInputFastq2 = null;


    /** Path to the genome digest file produced by GOPHER.*/
    @Parameter(names={"-d","--digest-file"}, required = true,description = "path to GOPHER digest file")
    private String digestFile;

    //** if this is set, an extra BAM file containg the rejected read pairs will be created */
    @Parameter(names={"-j", "--bad"}, description = "output bad (rejected) reads to separated file")
    private boolean outputRejectedReads;
    @Parameter(names={"-p", "--thread-num"},description = "number of threads used by bowtie2")
    private int threadNum = 1;
    @Parameter(names={"-l", "--lower-frag-size-limit"},required = true,description = "lower limit for fragment size")
    private int lowerFragSize = 50;
    @Parameter(names={"-u", "--upper-frag-size-limit"},required = true,description = "upper limit for fragment size")
    private int upperFragSize = 800;
    @Parameter(names={"-s", "--self-ligation-frag-size-limit"},required = true,description = "upper limit for self-ligation fragment size")
    private int upperSelfLigationFragSize = 3000;
    @Parameter(names={"-bsu","--bowtie-stringent-unique"}, description = "use stringent settings for definition of uniquely mapped reads")
    private boolean useStringentUniqueSettings = false;



    public AlignCommand(){}



    public void execute() throws DiachromaticException {

        String samFile1 = String.format("%s_%s_1.sam", this.filenamePrefix, getRandomPrefix(7));
        String samFile2 = String.format("%s_%s_2.sam", this.filenamePrefix, getRandomPrefix(7));
        logger.trace(String.format("About to read digests from %s.",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);
        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq1,samFile1,this.threadNum);
            runner.run();
            Bowtie2Runner runner2 = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq2,samFile2,this.threadNum);
            runner2.run();
            Aligner pairer = new Aligner(samFile1,samFile2, outputRejectedReads,filenamePrefix, digestMap, lowerFragSize, upperFragSize, upperSelfLigationFragSize, filenamePrefix,useStringentUniqueSettings);
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


    public String getRandomPrefix(int len) {
        char[] ch = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
                'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
                'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                'w', 'x', 'y', 'z' };

        char[] c=new char[len];
        Random random=new Random();
        for (int i = 0; i < len; i++) {
            c[i]=ch[random.nextInt(ch.length)];
        }
        return new String(c);
    }
}