package org.jax.diachromatic;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.google.common.collect.ImmutableSet;
import org.jax.diachromatic.command.*;
import org.jax.diachromatic.exception.DiachromaticException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * An application to process Hi-C data for differential reads counts in fragments surrounding the
 * transcription start site using probe design by VPV.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.1 (2017-11-15)
 */
public class Diachromatic {
    private static final Logger logger = LogManager.getLogger();

    @Parameter(names = {"-h", "--help"}, help = true, arity = 0, description = "display this help message")
    private boolean usageHelpRequested;

    private static final ImmutableSet<String> commandnames = ImmutableSet.of("truncate", "align", "count", "summarize");


    public static void main(String[] args) throws DiachromaticException {

        Diachromatic diachromatic = new Diachromatic();
        AlignCommand align = new AlignCommand();
        CountCommand count = new CountCommand();
        TruncateCommand truncate = new TruncateCommand();
        SummarizeCommand summarize = new SummarizeCommand();

        JCommander jc = JCommander.newBuilder().addObject(diachromatic).addCommand("truncate", truncate)
                .addCommand("align", align).addCommand("count", count).addCommand("summarize", summarize)
                .build();
        jc.setProgramName("java -jar Diachromatic.jar");
        try {
            jc.parse(args);
        } catch (ParameterException e) {
            // Note that by default, JCommand is OK with -h download but
            // not with download -h
            // The following hack makes things work with either option.
            String mycommand = null;
            String commandstring = String.join(" ", args);
            for (String a : args) {
                if (commandnames.contains(a)) {
                    mycommand = a;
                }
                if (a.equals("h")) {
                    if (mycommand != null) {
                        jc.usage(mycommand);
                    } else {
                        jc.usage();
                    }
                    System.exit(1);
                }
            }
            if (commandstring == null) { // user ran without any command
                jc.usage();
                System.exit(0);
            }
            logger.error("[ERROR] " + e.getMessage());
            logger.error("[ERROR] your command: " + commandstring);
            logger.error("[ERROR] enter java -jar Lr2pg -h for more information.");
            System.exit(1);
        }
        String parsedCommand = jc.getParsedCommand();
        if (parsedCommand == null) {
            jc.usage(); // user ran program with no arguments, probably help is want is wanted.
            System.exit(0);
        }
        if (!commandnames.contains(parsedCommand)) {
            logger.error("[ERROR] did not recognize command \"" + parsedCommand + "\"");
            logger.error("[ERROR] available commands are " + String.join(", ", commandnames));
            logger.error("[ERROR] enter java -jar Lr2pg -h for more information.");
            System.exit(1);
        }

        if (diachromatic.usageHelpRequested) {
            if (parsedCommand == null) {
                jc.usage();
            } else {
                jc.usage(parsedCommand);
            }
            System.exit(1);
        }

        if (jc.getParsedCommand() == null) {
            logger.error("[ERROR] no command passed");
            jc.usage();
            System.exit(1);
        }

        Command diachromaticCommand;
        switch (parsedCommand) {
            case "truncate":
                diachromaticCommand = truncate;
                break;
            case "align":
                diachromaticCommand = align;
                break;
            case "count":
                diachromaticCommand = count;
                break;
            case "summarize":
                diachromaticCommand = summarize;
                break;
            default:
                throw new DiachromaticException("Did not recognize command: " + parsedCommand);
        }
        try {
            diachromaticCommand.execute();
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }

    }

    public static String getVersion() {
        String version = "0.0.0";// default, should be overwritten by the following.
        try {
            Package p = Diachromatic.class.getPackage();
            version = p.getImplementationVersion();
        } catch (Exception e) {
            // do nothing
        }
        return version;
    }
}
