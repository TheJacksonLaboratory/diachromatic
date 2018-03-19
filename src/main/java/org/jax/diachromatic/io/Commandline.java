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

/**
 * Class to capture options and command from the command line.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.2 (2018-01-05)
 */
public class Commandline {

    private Command command=null;
    /** The default name of the file that is produced by the {@code digest} command. */
    private final static String DEFAULT_DIGEST_FILE_NAME="diachromaticDigest.txt";

    private final static String DEFAULT_OUTPUT_DIRECTORY="results";

    private final static String DEFAULT_TRUNCATION_SUFFIX="truncated";

    private final static String DEFAULT_OUTPUT_BAM_NAME="diachromatic-processed";
    /** Default size of margin of fragments used for calculating GC and repeat content. */
    private final static int DEFAULT_MARGIN_SIZE=250;

    /** Absolute path to the combined genome fasta file (which will be indexed only if necessary). */
    private String genomeFastaFile=null;
    private String outputFilePath=null;
    private String outputDirectory=null;
    private String enzyme=null;
    /** path to the bowtie2 executable. */
    private String bowtiepath=null;
    private String pathToBowtieIndex=null;
    private String pathToInputFastq1 =null;
    private String pathToInputFastq2 =null;
    private String pathToDiachromaticDigestFile=null;
    private String suffix=null;
    private boolean outputRejectedReads=false;
    private int marginsize;

    public Commandline(String args[]) {
        final CommandLineParser cmdLineGnuParser = new DefaultParser();

        final Options gnuOptions = constructOptions();
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
            if (category.length ==0) {
                printUsage("command missing");
            } else if (category.length > 1) {
                String cmd= Arrays.stream(category).collect(Collectors.joining("\t"));
                System.out.println(String.format("%d arguments for command: %s",category.length,cmd ));
            } else {
                mycommand=category[0];

            }
            if (commandLine.getArgs().length<1) {
                printUsage("no arguments passed");
                return;
            }
            if (commandLine.hasOption("b")) {
                this.bowtiepath=commandLine.getOptionValue("b");
            }
            if (commandLine.hasOption("d")) {
                this.pathToDiachromaticDigestFile=commandLine.getOptionValue("d");
            }
            if (commandLine.hasOption("e")) {
                this.enzyme=commandLine.getOptionValue("e");
            }
            if (commandLine.hasOption("g")) {
                this.genomeFastaFile=commandLine.getOptionValue("g");
            }
            if (commandLine.hasOption("i")) {
                this.pathToBowtieIndex=commandLine.getOptionValue("i");
            }
            if (commandLine.hasOption("j")) {
                outputRejectedReads=true;
            } else {
                outputRejectedReads=false;
            }
            if (commandLine.hasOption("m")) {
                String m = commandLine.getOptionValue("m");
                try {
                    this.marginsize = Integer.parseInt(m);
                } catch (NumberFormatException e) {
                    printUsage("[ERROR] -m option requires integer value. You passed \""+m+"\"");
                }
            } else {
                this.marginsize=DEFAULT_MARGIN_SIZE;
            }
            if (commandLine.hasOption("o")) {
                this.outputFilePath=commandLine.getOptionValue("o");
            }
            if (commandLine.hasOption("outdir")) {
                outputDirectory=commandLine.getOptionValue("outdir");
            }
            if (commandLine.hasOption("q")) {
                this.pathToInputFastq1 =commandLine.getOptionValue("q");
            }
            if (commandLine.hasOption("r")) {
                this.pathToInputFastq2 =commandLine.getOptionValue("r");
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
                if (this.genomeFastaFile == null) {
                    printUsage("-g option required for digest command");
                }
                if (this.enzyme==null) {
                    printUsage("-e option required for digest command");
                }
                if (this.outputFilePath == null) {
                    outputFilePath=DEFAULT_DIGEST_FILE_NAME;
                }
                this.command = new DigestCommand(this.genomeFastaFile, enzyme,this.outputFilePath,this.marginsize);

            } else if (mycommand.equalsIgnoreCase("truncate")) {
                if (this.outputDirectory == null) {
                    this.outputDirectory=DEFAULT_OUTPUT_DIRECTORY;
                } else if (this.pathToInputFastq1 == null) {
                    printUsage("-q option required for truncate command");
                } else if (this.pathToInputFastq2 == null) {
                    printUsage("-r option required for truncate command");
                } else if (this.enzyme == null) {
                    printUsage("-e option required for truncate command");
                }
                if (suffix==null) {
                    suffix=DEFAULT_TRUNCATION_SUFFIX;
                }
                //String outdir, String file1, String file2, String enzymeName
                this.command = new TruncateCommand(outputDirectory, pathToInputFastq1, pathToInputFastq2, enzyme,suffix);
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
                this.command=new MapCommand(bowtiepath,
                        pathToBowtieIndex,
                        pathToInputFastq1,
                        pathToInputFastq2,
                        outputFilePath,
                        pathToDiachromaticDigestFile,
                        outputRejectedReads);
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
    private static Options constructOptions()
    {
        final Options options = new Options();
        options.addOption("b", "bowtie", true, "path to bowtie2")
                .addOption("d", "digest", true, "path to diachromatic digest file")
                .addOption("e", "enzyme", true, "restriction enzyme name")
                .addOption("g", "genome", true, "path to genome FASTA file (with all chromosomes)")
                .addOption("i", "bowtieindex", true, "path to bowtie2 index")
                .addOption("j", "bad", false, "output bad (reJected) reads to separated file")
                .addOption("m","margin", true,"margin size for calculating GC and repeat content (default: 250 bp)")
                .addOption("o", "out", true, "name/path of output file/directory")
                .addOption("q", "q", true, "path to forward FASTQ input file")
                .addOption("r", "r", true, "path to reverse FASTQ input file")
                .addOption("s", "suffix", true, "suffix for output filenames")
                .addOption("outdir", "outdir", true, "path to output directory")
                .addOption( Option.builder( "f1" ).longOpt("file1").desc("path to fastq file 1").hasArg(true).argName("file1").build())
                .addOption( Option.builder( "f2" ).longOpt("file2").desc("path to fastq file 2").hasArg(true).argName("file2").build());
        return options;
    }

    public static String getVersion() {
        String version="0.0.0";// default, should be overwritten by the following.
        try {
            Package p = Commandline.class.getPackage();
            version = p.getImplementationVersion();
        } catch (Exception e) {
            // do nothing
        }
        return version;
    }

    /**
     * Print usage information to provided OutputStream.
     */
    private static void printUsage(String message)
    {

        String version=getVersion();
        final PrintWriter writer = new PrintWriter(System.out);
        writer.println(message);
        writer.println();
        writer.println("Program: Diachromatic (Analysis of Differential Capture Hi-C Interactions)");
        writer.println("Version: "+version);
        writer.println();
        writer.println("Usage: java -jar Diachromatic.jar <command> [options]");
        writer.println();
        writer.println("Available commands:");
        writer.println();
        writer.println("digest:");
        writer.println("\tjava -jar Diachromatic.jar digest -g <path> -e <enzyme> [-m <margin>] [-o <outfile>]");
        writer.println("\t<path>: path to a directory containing indexed genome FASTA files");
        writer.println("\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)");
        writer.println("\t<margin>: margin size in basepairs (Default: 250)");
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s\")",DEFAULT_DIGEST_FILE_NAME));
        writer.println();
        writer.println("truncate:");
        writer.println("\tjava -jar Diachromatic.jar truncate -q <forward.fq.gz> \\");
        writer.println("\t\t-r <reverse.fq.gz> -e <enzyme> -s <suffix> --outdir <directory>");
        writer.println("\t<forward.fq.gz>: path to the forward FASTQ file (may or may not be compressed with gzip)");
        writer.println("\t<reverse.fq.gz>: path to the reverse FASTQ file (may or may not be compressed with gzip)");
        writer.println("\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)");
        writer.println("\t<suffix>: suffix that will be added to the output truncated FASTQ files");
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s\")",DEFAULT_TRUNCATION_SUFFIX));
        writer.println();
        writer.println("map:");
        writer.println("\tjava -jar Diachromatic.jar map -b <bowtie2> -i <bowtie2-index> \\");
        writer.println("\t\t-q <forward.truncated.fq.gz> -r <reverse.truncated.fq.gz> \\");
        writer.println("\t\t-d <digest> [-o <outfile>] [-b]");
        writer.println("\t<bowtie2>: path to bowtie2 executable");
        writer.println("\t<bowtie2-index>: path to bowtie2 index for digested genome");
        writer.println("\t<forward.truncated.fq.gz>: path to the truncated forward FASTQ file");
        writer.println("\t<reverse.truncated.fq.gz>: path to the truncated reverse FASTQ file");
        writer.println("\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)");
        writer.println("\t<digest>: path to the digest file produced by the digest command");
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s.bam\")",DEFAULT_OUTPUT_BAM_NAME));
        writer.println("\t-b: output rejected reads to file (false if no -b option passed)");
        writer.println();
        writer.close();
        System.exit(0);
    }

}
