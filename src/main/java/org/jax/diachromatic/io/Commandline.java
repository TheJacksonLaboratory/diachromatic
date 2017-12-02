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

    private final static String DEFAULT_DIGEST_FILE_NAME="hicupCloneDigest.txt";

    private final static String DEFAULT_TRUNCATION_SUFFIX="truncated";



    private final static String DEFAULT_OUTPUT_BAM_NAME="diachromatic-processed";


    private String genomeDirectory=null;
    private String outputFilePath=null;
    private String file1=null;
    private String file2=null;
    private String enzyme=null;
    private String bowtiepath=null;
    private String pathToBowtieIndex=null;
    private String pathToInputFastq1 =null;
    private String pathToInputFastq2 =null;
    private String pathToDiachromaticDigestFile=null;
    private String suffix=null;

    public Commandline(String args[]) {
        final CommandLineParser cmdLineGnuParser = new DefaultParser();

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
            if (commandLine.hasOption("d")) {
                this.pathToDiachromaticDigestFile=commandLine.getOptionValue("d");
            }
            if (commandLine.hasOption("i")) {
                this.pathToBowtieIndex=commandLine.getOptionValue("i");
            }
            if (commandLine.hasOption("q")) {
                this.pathToInputFastq1 =commandLine.getOptionValue("q");
            }
            if (commandLine.hasOption("q2")) {
                this.pathToInputFastq2 =commandLine.getOptionValue("q2");
            }
            if (commandLine.hasOption("file1")) {
                this.file1 =commandLine.getOptionValue("file1");
            }
            if (commandLine.hasOption("file2")) {
                this.file2=commandLine.getOptionValue("file2");
            }
            if (commandLine.hasOption("s")) {
                this.suffix=commandLine.getOptionValue("s");
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
                if (this.enzyme==null) {
                    printUsage("-e option required for digest command");
                }
                if (this.outputFilePath == null) {
                    outputFilePath=DEFAULT_DIGEST_FILE_NAME;
                }
                this.command = new DigestCommand(this.genomeDirectory, enzyme,this.outputFilePath);

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
                if (suffix==null) {
                    suffix=DEFAULT_TRUNCATION_SUFFIX;
                }
                //String outdir, String file1, String file2, String enzymeName
                this.command = new TruncateCommand(outputFilePath, file1, file2, enzyme,suffix);
            } else if (mycommand.equalsIgnoreCase("map")) {
                if (this.bowtiepath==null) {
                    printUsage("-b option required for map command");
                }
                if (this.pathToBowtieIndex==null) {
                    printUsage("-i option (bowtie index) required for map command");
                }
                if (this.pathToInputFastq1 ==null) {
                    printUsage("-q option (FASTQ 1) required for map command");
                }
                if (this.pathToInputFastq2 ==null) {
                    printUsage("-r option (FASTQ 2) required for map command");
                }
                if (this.outputFilePath==null) {
                    outputFilePath=DEFAULT_OUTPUT_BAM_NAME;
                }
                if (pathToDiachromaticDigestFile==null) {
                    printUsage("-d option required for map command");
                }
                this.command=new MapCommand(bowtiepath,pathToBowtieIndex, pathToInputFastq1,pathToInputFastq2,outputFilePath,pathToDiachromaticDigestFile);
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
                .addOption("d", "digest", true, "path to diachromatic digest file")
                .addOption("s", "suffix", true, "suffix for output filenames")
                .addOption("i", "bowtieindex", true, "path to bowtie2 index")
                .addOption("q", "q1", true, "path to FASTQ input file")
                .addOption("q2", "q2", true, "path to FASTQ input file")
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
