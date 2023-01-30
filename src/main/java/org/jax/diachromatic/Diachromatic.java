package org.jax.diachromatic;




import org.jax.diachromatic.command.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;


import java.util.concurrent.Callable;

/**
 * An application to process Hi-C data for differential reads counts in fragments surrounding the
 * transcription start site using probe design by VPV.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.1 (2017-11-15)
 */
@CommandLine.Command(name = "diachromatic", mixinStandardHelpOptions = true, version = "0.2.9",
        description = "Hi-C and Capture Hi-C analysis")
public class Diachromatic implements Callable<Integer>  {
    private static final Logger logger = LoggerFactory.getLogger(Diachromatic.class);

    public static void main(String[] args) {
        if (args.length == 0) {
            // if the user doesn't pass any command or option, add -h to show help
            args = new String[]{"-h"};
        }
        long startTime = System.currentTimeMillis();
        CommandLine cline = new CommandLine(new Diachromatic())
                .addSubcommand("truncate", new TruncateCommand())
                .addSubcommand("align", new AlignCommand())
                .addSubcommand("count", new CountCommand())
                .addSubcommand("summarize", new SummarizeCommand());
        cline.setToggleBooleanFlags(false);
        int exitCode = cline.execute(args);
        long stopTime = System.currentTimeMillis();
        int elapsedTime = (int)((stopTime - startTime)*(1.0)/1000);
        if (elapsedTime > 3599) {
            int elapsedSeconds = elapsedTime % 60;
            int elapsedMinutes = (elapsedTime/60) % 60;
            int elapsedHours = elapsedTime/3600;
            logger.info("Elapsed time {}h:{}m-{}m",elapsedHours,elapsedMinutes,elapsedSeconds);
        }
        else if (elapsedTime>59) {
            int elapsedSeconds = elapsedTime % 60;
            int elapsedMinutes = (elapsedTime/60) % 60;
            logger.info("Elapsed time {} min, {} sec",elapsedMinutes,elapsedSeconds);
        } else {
            logger.info("Elapsed time {} seconds.", (stopTime - startTime) * (1.0) / 1000 );
        }
        System.exit(exitCode);
    }

    @Override
    public Integer call() {
        // work done in subcommands
        return 0;
    }

    public static String getVersion() {
        String version = "0.7.0";// default, should be overwritten by the following.
        try {
            Package p = Diachromatic.class.getPackage();
            version = p.getImplementationVersion();
        } catch (Exception e) {
            // do nothing
        }
        return version;
    }
}
