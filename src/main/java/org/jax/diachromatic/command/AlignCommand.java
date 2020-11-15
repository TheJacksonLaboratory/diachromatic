package org.jax.diachromatic.command;




import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import org.jax.diachromatic.align.DigestMap;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.align.Aligner;
import org.jax.diachromatic.align.Bowtie2Runner;
import picocli.CommandLine;

import java.io.File;
import java.io.IOException;
import java.util.Random;
import java.util.concurrent.Callable;

/**
 * Class to coordinate bowtie2-mapping of the truncated FASTQ files followed by Q/C and filtering of the mapped reads.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.3 (2018-03-20)
 */
@CommandLine.Command(name = "align",
        aliases = {"A"},
        mixinStandardHelpOptions = true,
        description = "align with uses bowtie2 to align Hi-C reads and then performs Q/C, artifact filtering and pairing of valid read pairs.")
public class AlignCommand extends Command implements Callable<Integer>  {
    private static final Logger logger = LoggerFactory.getLogger(AlignCommand.class);
    /** Path to the bowtie2 executable, e.g., {@code /usr/bin/bowtie2}. */
    @CommandLine.Option(names={"-b","--bowtie-path"},required = true, description ="Path to bowtie2.", order = 1)
    private String bowtiepath;

    /** Path to the bowtie2 index files. Note that the index is made up of multiple files, e.g.,
     * hg19.1.bt2,  hg19.3.bt2,  hg19.rev.1.bt2, hg19.2.bt2,  hg19.4.bt2,  hg19.rev.2.bt2. Assuming all files
     * are in a directory called {@code /path/to/index/}, this parameter should be {@code /path/to/index/hg19}.*/
    @CommandLine.Option(names={"-i", "--bowtie-index"}, required = true, description ="Path to bowtie2 index.", order = 2)
    private String pathToBowtieIndex;

    @CommandLine.Option(names={"-p", "--thread-num"},description = "Number of threads used by bowtie2.", order = 3)
    private int threadNum = 1;

    @CommandLine.Option(names={"-bsu","--bowtie-stringent-unique"},
            description = "Use stringent settings for definition of uniquely mapped reads.", order = 4)
    private boolean useStringentUniqueSettings = false;

    /** Path to the forward truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    @CommandLine.Option(names={"-q","fastq-r1"}, required = true, description = "Path to truncated forward FASTQ input file.",order = 5)
    private String pathToInputFastq1 = null;

    /** Path to the reverse truncated FASTQ file produced by {@link org.jax.diachromatic.command.TruncateCommand}. */
    @CommandLine.Option(names={"-r","fastq-r2"}, required = true, description = "Path to truncated reverse FASTQ input file.", order = 6)
    private String pathToInputFastq2 = null;

    /** Path to the genome digest file produced by GOPHER.*/
    @CommandLine.Option(names={"-d","--digest-file"}, required = true, description = "Path to GOPHER digest file.", order = 7)
    private String digestFile;
    @CommandLine.Option(names={"-l", "--lower-frag-size-limit"}, description = "Lower limit for fragment size.", order = 8)
    private int lowerFragSize = 50;
    @CommandLine.Option(names={"-u", "--upper-frag-size-limit"}, description = "Upper limit for fragment size.", order = 9)
    private int upperFragSize = 800;
    @CommandLine.Option(names={"-s", "--self-ligation-frag-size-limit"}, description = "Upper limit for self-ligation fragment size.", order = 10)
    private int upperSelfLigationFragSize = 3000;
    /** if this is set, an extra BAM file containg the rejected read pairs will be created */
    @CommandLine.Option(names={"-j", "--bad"}, description = "Output bad (rejected) reads to separated file.", order = 11)
    private boolean outputRejectedReads=false;
    /** if this is set, an extra BAM file containg the rejected read pairs will be created */
    @CommandLine.Option(names={"-k", "--keep-sam"}, description = "Do not delete temporary SAM files.",order = 12)
    private boolean keepSamFiles=false;

    public AlignCommand(){}


    @Override
    public Integer call() throws DiachromaticException {

        makeOutdirectoryIfNeeded();

        String outputDirAndFilePrefix = String.format("%s%s%s", outputDir,File.separator,filenamePrefix);

        String samFile1 = String.format("%s_%s_1.sam", outputDirAndFilePrefix, getRandomPrefix(7));
        String samFile2 = String.format("%s_%s_2.sam", outputDirAndFilePrefix, getRandomPrefix(7));
        logger.trace(String.format("About to read digests from %s.",digestFile));
        DigestMap digestMap = new DigestMap(digestFile);
        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq1,samFile1,this.threadNum);
            runner.run();
            Bowtie2Runner runner2 = new Bowtie2Runner(bowtiepath,pathToBowtieIndex,pathToInputFastq2,samFile2,this.threadNum);
            runner2.run();

            Aligner pairer = new Aligner(samFile1,samFile2, outputRejectedReads, outputDirAndFilePrefix, digestMap, lowerFragSize, upperFragSize, upperSelfLigationFragSize, filenamePrefix,useStringentUniqueSettings);
            pairer.inputSAMfiles();
            pairer.printStatistics();
            if(!keepSamFiles) {
                File file = new File(samFile1);
                file.delete();
                file = new File(samFile2);
                file.delete();
            }
        } catch (DiachromaticException | IOException e){
            e.printStackTrace();
        }
        return 0;
    }
    @Override
    public String toString() {return "diachromatic:align";}


    /**
     * Needed for creation of random file names for temporary SAM files.
     *
     * @param len length of the desired random prefix
     * @return a random prefix that will be used for a temporary alignment file
     */
    private String getRandomPrefix(int len) {
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
