package org.jax.diachromatic.io;

import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.commons.cli.*;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.DigestCommand;
import org.jax.diachromatic.command.MapCommand;
import org.jax.diachromatic.command.TruncateCommand;
import org.jax.diachromatic.exception.DiachromaticException;


public class Commandline {

    private Command command=null;

    private String genomeDirectory=null;
    private String outputFilePath=null;
    private String digestfilename="hicupCloneDigest.txt";
    private String file1=null;
    private String file2=null;
    private String enzyme=null;
    private String bowtiepath=null;

    public Commandline(String args[]) {
        final CommandLineParser cmdLineGnuParser = new GnuParser();

        final Options gnuOptions = constructGnuOptions();
        org.apache.commons.cli.CommandLine commandLine;

        String mycommand = null;
        String clstring="";
        if (args!=null && args.length>0) {
            clstring= Arrays.stream(args).collect(Collectors.joining(" "));
        }
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
            if (commandLine.hasOption("e")) {
                this.enzyme=commandLine.getOptionValue("e");
            }
            if (commandLine.hasOption("b")) {
                this.bowtiepath=commandLine.getOptionValue("b");
            }
            if (commandLine.hasOption("file1")) {
                this.file1 =commandLine.getOptionValue("file1");
            }
            if (commandLine.hasOption("file2")) {
                this.file2=commandLine.getOptionValue("file2");
            }
        }
        catch (ParseException parseException)  // checked exception
        {
            String msg = String.format("Could not parse options %s [%s]",clstring, parseException.toString());
           printUsage(msg );
        }
        try {
            if (mycommand.equals("digest")) {
                if (this.genomeDirectory == null) {
                    printUsage("-g option required for digest command");
                }
                if (this.outputFilePath == null) {
                    this.command = new DigestCommand(this.genomeDirectory);
                } else {
                    this.command = new DigestCommand(this.genomeDirectory, this.outputFilePath);
                }
            } else if (mycommand.equalsIgnoreCase("truncate")) {
                if (this.outputFilePath == null) {
                    printUsage("-o option required for truncate command");
                } else if (this.file1 == null) {
                    printUsage("--file1 option required for truncate command");
                } else if (this.file2 == null) {
                    printUsage("--file2 option required for truncate command");
                } else if (this.enzyme == null) {
                    printUsage("-e option required for truncate command");
                }
                //String outdir, String file1, String file2, String enzymeName
                this.command = new TruncateCommand(outputFilePath, file1, file2, enzyme);
            } else if (mycommand.equalsIgnoreCase("map")) {
                if (this.bowtiepath==null) {
                    printUsage("-b option required for map command");
                }
                this.command=new MapCommand(bowtiepath);
            } else {
                printUsage(String.format("Did not recognize command: %s", mycommand));
            }
        } catch (DiachromaticException de) {
            de.printStackTrace();
            printUsage(de.getMessage());
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
        gnuOptions.addOption("o", "out", true, "name/path of output file/directory")
                .addOption("g", "genome", true, "genome directory (with FASTA files)")
                .addOption("e", "enzyme", true, "restriction enzyme name")
                .addOption("b", "bowtie", true, "path to bowtie2")
        .addOption( Option.builder( "f1" ).longOpt("file1").desc("path to fastq file 1").hasArg(true).argName("file1").build())
         .addOption( Option.builder( "f2" ).longOpt("file2").desc("path to fastq file 2").hasArg(true).argName("file2").build());
        return gnuOptions;
    }



    /**
     * Print usage information to provided OutputStream.
     */
    public static void printUsage(String message)
    {
        final PrintWriter writer = new PrintWriter(System.out);
        final HelpFormatter usageFormatter = new HelpFormatter();
        final String applicationName="java -jar diachromatic.jar command";
        final Options options=constructGnuOptions();
        writer.println(message);
        usageFormatter.printUsage(writer, 120, applicationName, options);
        writer.println("\twhere command is one of digest,truncate,....");
        writer.println("\t- digest -g genome [-o outputname]: Digest genome at directory (-g), output to file.");
        writer.println("\t- truncate -o outdir --file1 example1.fq.gz --file2 example2.fq.gz - enzymeName: truncate fastq files.");
        writer.close();
        System.exit(0);
    }

}
