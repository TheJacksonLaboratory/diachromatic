package org.jax.diachromatic;


import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import org.jax.diachromatic.command.AlignCommand;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.CountCommand;
import org.jax.diachromatic.command.TruncateCommand;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.io.Commandline;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
/**
 * An application to process Hi-C data for differential reads counts in fragments surrounding the
 * transcription start site using probe design by VPV.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.1 (2017-11-15)
 */
public class Diachromatic {
    private static final Logger logger = LogManager.getLogger();

    @Parameter(names = {"-h", "--help"}, help = true, arity = 0,description = "display this help message")
    private boolean usageHelpRequested;


    public static void main(String[] args) throws DiachromaticException {
//        Commandline clp = new Commandline(args);
//        Command command = clp.getCommand();
//        logger.trace(String.format("Running command %s",command));
//        command.execute();

        Diachromatic diachromatic = new Diachromatic();
        AlignCommand align = new AlignCommand();
        CountCommand count = new CountCommand();
        TruncateCommand truncate = new TruncateCommand();

        JCommander jc = JCommander.newBuilder()
                .addObject(diachromatic)
                .addCommand("truncate", truncate)
                .addCommand("align", align)
                .addCommand("count", count)
                .build();
        jc.setProgramName("java -jar diachromatic.jar");
        try {
            jc.parse(args);
        } catch (ParameterException e) {
            // Note that by default, JCommand is OK with -h download but
            // not with download -h
            // The following hack makes things work with either option.
            String commandString=null;
            for (String a:args) {
//                if (commandnames.contains(a)) {
//                    commandString=a;
//                }
                if (a.contains("h")) {
                    if (commandString!=null) {
                        jc.usage(commandString);
                    } else {
                        jc.usage();
                    }
                    System.exit(1);
                }
            }
            jc.usage();
            System.exit(1);
        }
        String parsedCommand = jc.getParsedCommand();

        if ( diachromatic.usageHelpRequested) {
            if (parsedCommand==null) {
                jc.usage();
            } else {
                jc.usage(parsedCommand);
            }
            System.exit(1);
        }

        if (jc.getParsedCommand()==null ) {
            System.err.println("[ERROR] no command passed");
            jc.usage();
            System.exit(1);
        }

        if ( diachromatic.usageHelpRequested) {

            jc.usage();
            System.exit(1);
        }
        String command = jc.getParsedCommand();
        Command diachromaticCommand=null;
        switch (command) {
            case "truncate":
                diachromaticCommand = truncate;
                break;
            case "align":
                diachromaticCommand = align;
                break;
            case "count":
                diachromaticCommand = count;
                break;
            default:
                throw new DiachromaticException("Did not recognize command: "+ command);
        }
        try {
            diachromaticCommand.execute();
        } catch (DiachromaticException e) {
            e.printStackTrace();
        }
    }
}
