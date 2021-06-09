package org.jax.diachromatic;




import org.jax.diachromatic.command.AlignCommand;
import org.jax.diachromatic.command.CountCommand;
import org.jax.diachromatic.command.SummarizeCommand;
import org.jax.diachromatic.command.TruncateCommand;
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
        description = "Structural variant annotation")
public class Diachromatic implements Callable<Integer>  {
    private static final Logger logger = LoggerFactory.getLogger(Diachromatic.class);


//    @Parameter(names = {"-h", "--help"}, help = true, arity = 0, description = "display this help message")
//    private boolean usageHelpRequested;
//
//    private static final ImmutableSet<String> commandnames = ImmutableSet.of("truncate", "align", "count", "summarize");


    public static void main(String[] args) {
        if (args.length == 0) {
            // if the user doesn't pass any command or option, add -h to show help
            args = new String[]{"-h"};
        }
        long startTime = System.currentTimeMillis();
        CommandLine cline = new CommandLine(new Diachromatic())
                .addSubcommand("align", new AlignCommand())
                .addSubcommand("count", new CountCommand())
                .addSubcommand("truncate", new TruncateCommand())
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


//        Diachromatic diachromatic = new Diachromatic();
//        AlignCommand align = new AlignCommand();
//        CountCommand count = new CountCommand();
//        TruncateCommand truncate = new TruncateCommand();
//        SummarizeCommand summarize = new SummarizeCommand();
//
//        JCommander jc = JCommander.newBuilder().addObject(diachromatic).addCommand("truncate", truncate)
//                .addCommand("align", align).addCommand("count", count).addCommand("summarize", summarize)
//                .build();
//        jc.setProgramName("java -jar Diachromatic.jar");
//        try {
//            jc.parse(args);
//        } catch (ParameterException e) {
//            // Note that by default, JCommand is OK with -h download but
//            // not with download -h
//            // The following hack makes things work with either option.
//            String mycommand = null;
//            String commandstring = String.join(" ", args);
//            for (String a : args) {
//                if (commandnames.contains(a)) {
//                    mycommand = a;
//                }
//                if (a.equals("h")) {
//                    if (mycommand != null) {
//                        jc.usage(mycommand);
//                    } else {
//                        jc.usage();
//                    }
//                    System.exit(1);
//                }
//            }
//            if (commandstring == null) { // user ran without any command
//                jc.usage();
//                System.exit(0);
//            }
//            logger.error("[ERROR] " + e.getMessage());
//            logger.error("[ERROR] your command: " + commandstring);
//            logger.error("[ERROR] enter java -jar Lr2pg -h for more information.");
//            System.exit(1);
//        }
//        String parsedCommand = jc.getParsedCommand();
//        if (parsedCommand == null) {
//            jc.usage(); // user ran program with no arguments, probably help is want is wanted.
//            System.exit(0);
//        }
//        if (!commandnames.contains(parsedCommand)) {
//            logger.error("[ERROR] did not recognize command \"" + parsedCommand + "\"");
//            logger.error("[ERROR] available commands are " + String.join(", ", commandnames));
//            logger.error("[ERROR] enter java -jar Lr2pg -h for more information.");
//            System.exit(1);
//        }
//
//        if (diachromatic.usageHelpRequested) {
//            if (parsedCommand == null) {
//                jc.usage();
//            } else {
//                jc.usage(parsedCommand);
//            }
//            System.exit(1);
//        }
//
//        if (jc.getParsedCommand() == null) {
//            logger.error("[ERROR] no command passed");
//            jc.usage();
//            System.exit(1);
//        }
//
//        Command diachromaticCommand;
//        switch (parsedCommand) {
//            case "truncate":
//                diachromaticCommand = truncate;
//                break;
//            case "align":
//                diachromaticCommand = align;
//                break;
//            case "count":
//                diachromaticCommand = count;
//                break;
//            case "summarize":
//                diachromaticCommand = summarize;
//                break;
//            default:
//                throw new DiachromaticException("Did not recognize command: " + parsedCommand);
//        }
//        try {
//            diachromaticCommand.execute();
//        } catch (DiachromaticException e) {
//            e.printStackTrace();
//        }

    }

    @Override
    public Integer call() {
        // work done in subcommands
        return 0;
    }

    public static String getVersion() {
        String version = "0.6.0";// default, should be overwritten by the following.
        try {
            Package p = Diachromatic.class.getPackage();
            version = p.getImplementationVersion();
        } catch (Exception e) {
            // do nothing
        }
        return version;
    }
}
