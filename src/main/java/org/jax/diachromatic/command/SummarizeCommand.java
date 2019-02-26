package org.jax.diachromatic.command;


import com.beust.jcommander.Parameter;
import com.beust.jcommander.Parameters;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.summarize.Summarizer;

/**
 * Writes an HTML page with a summary of the entire analysis. We expect to find the files with the
 * prefix as defined and used for {@link TruncateCommand}, {@link AlignCommand} and {@link CountCommand}.
 * TODO
 * define CLI
 * for now
 * java -jar Diachromatic.jar summarize -o tst -x test -t truncate_test/srr5633684.truncation.stats.txt
 */
@Parameters(commandDescription = "The summarize command outputs an HTML file with a summary of the analysis.")
public class SummarizeCommand extends Command {
    private static final Logger logger = LogManager.getLogger();

    /** Path to the genome digest file produced by GOPHER. */
    @Parameter(names={"-t","--truncate"}, required = true, description = "Path to diachromatic truncate summarize file.", order = 4)
    private String truncateFile = null;






    public SummarizeCommand() {
    }




    @Override
    public void execute() throws DiachromaticException {
        Summarizer summarizer = new Summarizer(truncateFile);
        summarizer.outputFile(this.filenamePrefix);
    }
}
