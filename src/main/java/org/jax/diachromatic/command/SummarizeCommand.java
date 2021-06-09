package org.jax.diachromatic.command;



import org.jax.diachromatic.summarize.Summarizer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;

import java.io.File;
import java.util.concurrent.Callable;

/**
 * Writes an HTML page with a summary of the entire analysis. We expect to find the files with the
 * prefix as defined and used for {@link TruncateCommand}, {@link AlignCommand} and {@link CountCommand}.
 */
@CommandLine.Command(name = "summarize",
        aliases = {"S"},
        mixinStandardHelpOptions = true,
        description = "The summarize command outputs an HTML file with a summary of the analysis.")
public class SummarizeCommand extends Command implements Callable<Integer> {
    private static final Logger logger = LoggerFactory.getLogger(SummarizeCommand.class);

    /** Path to text file with summary statistics produced with the truncate command. */
    @CommandLine.Option(names={"-t","--truncate"}, required = true, description = "Path to diachromatic truncate statistics file.", order = 2)
    private String truncateFile = null;

    /** Path to text file with summary statistics produced with the align command. */
    @CommandLine.Option(names={"-a","--align"}, required = true, description = "Path to diachromatic align statistics file.", order = 3)
    private String alignFile = null;

    /** Path to text file with summary statistics produced with the count command. */
    @CommandLine.Option(names={"-c","--count"}, required = true, description = "Path to diachromatic statistics output file.", order = 4)
    private String countFile = null;

    public SummarizeCommand() {
    }

    @Override
    public Integer call() {
        makeOutdirectoryIfNeeded();
        String outputDirAndFilePrefix = String.format("%s%s%s", outputDir, File.separator, filenamePrefix);
        Summarizer summarizer = new Summarizer(truncateFile, alignFile, countFile);
        summarizer.outputFile(outputDirAndFilePrefix);
        return 0;
    }
}
