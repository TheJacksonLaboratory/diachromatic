package org.jax.diachromatic.command;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.summarize.Summarizer;

import java.io.File;

/**
 * Writes an HTML page with a summary of the entire analysis. We expect to find the files with the
 * prefix as defined and used for {@link TruncateCommand}, {@link AlignCommand} and {@link CountCommand}.
 */
@Parameters(commandDescription = "The summarize command outputs an HTML file with a summary of the analysis.")
public class SummarizeCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to text file with summary statistics produced with the truncate command. */
    @Parameter(names={"-t","--truncate"}, required = true, description = "Path to diachromatic truncate statistics file.", order = 2)
    private String truncateFile = null;

    /** Path to text file with summary statistics produced with the align command. */
    @Parameter(names={"-a","--align"}, required = true, description = "Path to diachromatic align statistics file.", order = 3)
    private String alignFile = null;

    /** Path to text file with summary statistics produced with the count command. */
    @Parameter(names={"-c","--count"}, required = true, description = "Path to diachromatic statistics output file.", order = 4)
    private String countFile = null;

    public SummarizeCommand() {
    }

    @Override
    public void execute() throws DiachromaticException {
        makeOutdirectoryIfNeeded();
        String outputDirAndFilePrefix = String.format("%s%s%s", outputDir, File.separator, filenamePrefix);
        Summarizer summarizer = new Summarizer(truncateFile, alignFile, countFile);
        summarizer.outputFile(outputDirAndFilePrefix);
    }
}
