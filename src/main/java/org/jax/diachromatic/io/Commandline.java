package org.jax.diachromatic.io;

import java.io.PrintWriter;

import org.apache.commons.cli.*;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.DigestCommand;


public class Commandline {

    private Command command=null;

    private String genomeDirectory=null;
    private String outputFilePath=null;

    public Commandline(String args[]) {
        final CommandLineParser cmdLineGnuParser = new GnuParser();

        final Options gnuOptions = constructGnuOptions();
        org.apache.commons.cli.CommandLine commandLine;

        String mycommand = null;
        try
        {
            commandLine = cmdLineGnuParser.parse(gnuOptions, args);
            String category[] = commandLine.getArgs();
            if (category.length != 1) {
                printUsage("command missing");
            } else {
                mycommand=category[0];

            }
            if (commandLine.getArgs().length<1) {
                printUsage("no arguments passed");
                return;
            }
            if (commandLine.hasOption("g")) {
                this.genomeDirectory=commandLine.getOptionValue("g");
            }
            if (commandLine.hasOption("o")) {
                this.outputFilePath=commandLine.getOptionValue("o");
            }
        }
        catch (ParseException parseException)  // checked exception
        {
            System.err.println(
                    "Encountered exception while parsing using GnuParser:\n"
                            + parseException.getMessage() );
        }
        if (mycommand.equals("digest")) {
            if (this.genomeDirectory==null) {
                printUsage("-g option required for digest command");
            }
            if (this.outputFilePath==null) {
                this.command = new DigestCommand(this.genomeDirectory);
            } else {
                this.command = new DigestCommand(this.genomeDirectory,this.outputFilePath);
            }
        } else {
            printUsage(String.format("Did not recognize command: %s", mycommand));
        }
    }


    public Command getCommand() {
        return command;
    }

    /**
     * Construct and provide GNU-compatible Options.
     *
     * @return Options expected from command-line of GNU form.
     */
    public static Options constructGnuOptions()
    {
        final Options gnuOptions = new Options();
        gnuOptions.addOption("o", "out", true, "name/path of output file")
                .addOption("g", "genome", true, "genome directory (with FASTA files)");
                /*.addOption("n", true, "Number of copies")
                .addOption("c","cpe",true,"Symbol of core promoter element to be analyzed")
                .addOption( OptionBuilder.withLongOpt( "cpe1" )
                        .withDescription( "core promoter element 1" )
                        .hasArg()
                        .create())
                .addOption( OptionBuilder.withLongOpt( "cpe2" )
                        .withDescription( "core promoter element 2" )
                        .hasArg()
                        .create());*/
        return gnuOptions;
    }



    /**
     * Print usage information to provided OutputStream.
     */
    public static void printUsage(String message)
    {
        final PrintWriter writer = new PrintWriter(System.out);
        final HelpFormatter usageFormatter = new HelpFormatter();
        final String applicationName="Oncembobulator";
        final Options options=constructGnuOptions();
        usageFormatter.printUsage(writer, 80, applicationName, options);
        writer.print("\t where command is one of download-oncokb, download-clinvar, parse.\n");
        writer.print("\t download-oncokb: Download the OnkoKB data to the data directory.\n");
        writer.print("\t download-clinvar: Download the ClinVar data to the data directory.\n");
        writer.print("\t jannovar: Download hg38 transcript data and create Jannovar transcript file.\n");
        writer.print("\t undiscombobulate: map the the mutations.\n");
        writer.close();
        System.exit(0);
    }

}
