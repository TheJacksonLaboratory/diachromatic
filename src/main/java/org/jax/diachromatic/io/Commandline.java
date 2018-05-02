package org.jax.diachromatic.io;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

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
    private static final Logger logger = LogManager.getLogger();

    private Command command=null;
    /** The default name of the file that is produced by the {@code digest} command. */
    private final static String DEFAULT_DIGEST_FILE_NAME="default";

    private final static String DEFAULT_OUTPUT_DIRECTORY="results";
    private final static String DEFAULT_FILENAME_PREFIX="prefix";


    private final static String DEFAULT_TRUNCATION_SUFFIX="truncated";
    private final static String DEFAULT_MAPPING_SUFFIX="mapped";
    



    private final static String DEFAULT_OUTPUT_BAM_NAME="diachromatic-processed";
    /** Default size of margin of fragments used for calculating GC and repeat content. */
    private final static int DEFAULT_MARGIN_SIZE=250;

    /** Absolute path to the combined genome fasta file (which will be indexed only if necessary). */
    private String genomeFastaFile=null;
    private String outputFilePath=null;
    private String outputDirectory=DEFAULT_OUTPUT_DIRECTORY;
    private String enzyme=null;
    /** path to the bowtie2 executable. */
    private String bowtiepath=null;
    private String pathToBowtieIndex=null;
    private String pathToInputFastq1 =null;
    private String pathToInputFastq2 =null;
    private String pathToDiachromaticDigestFile=null;
    private String pathToActiveDigestsFile=null;
    private String outprefix=null;
    private boolean outputRejectedReads=false;
    private boolean doHelp=false;
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
            // parse commandline
            commandLine = cmdLineGnuParser.parse(gnuOptions, args);
            String category[] = commandLine.getArgs();
            if (category.length == 0) {
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
            if (commandLine.hasOption("a")) {
                this.pathToActiveDigestsFile=commandLine.getOptionValue("a");
            }
            if (commandLine.hasOption("e")) {
                this.enzyme=commandLine.getOptionValue("e");
            }
            if (commandLine.hasOption("g")) {
                this.genomeFastaFile=commandLine.getOptionValue("g");
            }
            if (commandLine.hasOption("h")) {
                this.doHelp=true;
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
            if (commandLine.hasOption("outprefix")) {
                this.outprefix=commandLine.getOptionValue("outprefix");
            }

            // create directory for output in any case
            makeOutdirectoryIfNeeded();
        }
        catch (ParseException parseException)  // checked exception
        {
            String msg = String.format("Could not parse options %s [%s]",clstring, parseException.toString());
           printUsage(msg );
        }
        if (doHelp ) {
            switch (mycommand) {
                case "digest": printHelpHeader(); printDigestHelp(true); break;
                case "truncate": printHelpHeader(); printTruncateHelp(true); break;
                case "map": printHelpHeader(); printMapHelp(true); break;
            }
            System.exit(1);
        }
        try {
            if (mycommand.equals("digest")) {
                if (this.genomeFastaFile == null) {
                    printDigestHelp("-g option required for digest command");
                }
                if (this.enzyme==null) {
                    printDigestHelp("-e option required for digest command");
                }
                if (this.outputFilePath == null) {
                    outputFilePath=DEFAULT_DIGEST_FILE_NAME;
                }
                this.command = new DigestCommand(this.genomeFastaFile, enzyme,this.outputFilePath,this.marginsize);

            } else if (mycommand.equalsIgnoreCase("truncate")) {
                if (this.outputDirectory == null) {
                    this.outputDirectory=DEFAULT_OUTPUT_DIRECTORY;
                } else if (this.pathToInputFastq1 == null) {
                    printTruncateHelp("-q option required for truncate command");
                } else if (this.pathToInputFastq2 == null) {
                    printTruncateHelp("-r option required for truncate command");
                } else if (this.enzyme == null) {
                    printTruncateHelp("-e option required for truncate command");
                }
                if (outprefix==null) {
                    outprefix=DEFAULT_FILENAME_PREFIX;
                }
                //String outdir, String file1, String file2, String enzymeName
                this.command = new TruncateCommand(pathToInputFastq1, pathToInputFastq2, enzyme, outputDirectory, outprefix);
            } else if (mycommand.equalsIgnoreCase("map")) {
                if (this.bowtiepath==null) {
                    printMapHelp("-b option required for map command");
                }
                if (this.pathToBowtieIndex==null) {
                    printMapHelp("-i option (bowtie index) required for map command");
                }
                if (this.pathToInputFastq1 ==null) {
                    printMapHelp("-q option (FASTQ 1) required for map command");
                }
                if (this.pathToInputFastq2 ==null) {
                    printMapHelp("-r option (FASTQ 2) required for map command");
                }
                if (this.outputFilePath==null) {
                    outputFilePath=DEFAULT_OUTPUT_BAM_NAME;
                }
                if (pathToDiachromaticDigestFile==null) {
                    printMapHelp("-d option required for map command");
                }
                this.command=new MapCommand(
                        bowtiepath,
                        pathToBowtieIndex,
                        pathToInputFastq1,
                        pathToInputFastq2,
                        outputFilePath,
                        pathToDiachromaticDigestFile,
                        pathToActiveDigestsFile,
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
                .addOption("a", "active-digests", true, "path to BED file with active digests")
                .addOption("e", "enzyme", true, "restriction enzyme name")
                .addOption("g", "genome", true, "path to genome FASTA file (with all chromosomes)")
                .addOption("h", "help", false, "shows help for current command")
                .addOption("i", "bowtieindex", true, "path to bowtie2 index")
                .addOption("j", "bad", false, "output bad (reJected) reads to separated file")
                .addOption("m","margin", true,"margin size for calculating GC and repeat content (default: 250 bp)")
                .addOption("o", "out", true, "name/path of output file/directory")
                .addOption("q", "q", true, "path to forward FASTQ input file")
                .addOption("r", "r", true, "path to reverse FASTQ input file")
                .addOption("outdir", "outdir", true, "path to output directory")
                .addOption("outprefix", "outprefix", true, "outprefix for files in output directory")
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


    private static void printDigestHelp(String message) {
        System.out.println("\n"+ message + "\n");
        printDigestHelp(true);
        System.out.println();
        System.exit(0);
    }

    private static void printDigestHelp(boolean specific) {

        if (specific) {
            System.out.println("Command and options for digest function");
            System.out.println("\tdigest creates an in silico digest of the genome that is need in later steps of the analysis");
        }

        System.out.println("digest:\n" +
            "\tjava -jar Diachromatic.jar digest -g <path> -e <enzyme> [-o <outfile>] [-m <margin>]\n\n" +

            "\t\t<path>: path to a directory containing indexed genome FASTA files (e.g. genomes/hg19.fa)\n" +
            "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)\n" +
            "\t\t<margin>: margin size in basepairs ? (Default: 250)\n" +
            "\t\t<outfile>: optional name of output file (Default: <genome>_<enzyme>_digest.tsv)\n\n");
    }

    private static void printTruncateHelp(String message) {
        System.out.println("\n"+ message + "\n");
        printTruncateHelp(true);
        System.out.println();
        System.exit(0);
    }

    private static void printTruncateHelp(boolean specific) {
        if (specific) {
            System.out.println("Command and options for truncate function");
            System.out.println("\ttruncate searches for Hi-C religation sequences and truncates reads accordingly");
        }
        System.out.println("truncate:\n" +
             "\tjava -jar Diachromatic.jar truncate -q <forward.fq.gz> \\ \n"+
            "\t\t\t-r <reverse.fq.gz> -e <enzyme> -outdir <directory> -outprefix <filename_prefix>\n\n"+
                "\t\t<forward.fq.gz>: path to the forward FASTQ file (may or may not be compressed with gzip)\n"+
        "\t\t<reverse.fq.gz>: path to the reverse FASTQ file (may or may not be compressed with gzip)\n"+
        "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)\n"+
        "\t\t<directory>: directory containing the output of the truncate command (Default: " + DEFAULT_OUTPUT_DIRECTORY + ")\n"+
        "\t\t<filename_prefix>: prefix for all generated files in output directory (Default: " + DEFAULT_FILENAME_PREFIX + ")\n");
        //String.format("\t\t<outfile>: optional name of output file (Default: \"%s\")",DEFAULT_TRUNCATION_SUFFIX));
    }

    private static void printMapHelp(String message) {
        System.out.println("\n"+ message + "\n");
        printMapHelp(true);
        System.out.println();
        System.exit(0);
    }


    private static void printMapHelp(boolean specific) {
        if (specific) {
            System.out.println("Command and options for map function");
            System.out.println("\tmap map uses bowtie2 to map reads and then performs Q/C and repairing.");
        }
        System.out.println("map:\n" +
        "\tjava -jar Diachromatic.jar map -b <bowtie2> -i <bowtie2-index> \\ \n" +
        "\t\t\t-q <forward.truncated.fq.gz> -r <reverse.truncated.fq.gz> \\ \n" +
        "\t\t\t-d <digest> [-o <outfile>] [-b <output-rejected>]\n\n" +
        "\t\t\t[-a <active-digests>] [-o <outfile>] [-b]\n" +
        "\t\t<bowtie2>: path to bowtie2 executable\n" +
        "\t\t<bowtie2-index>: path to bowtie2 index for digested genome\n" +
        "\t\t<forward.truncated.fq.gz>: path to the truncated forward gzipped FASTQ file\n" +
        "\t\t<reverse.truncated.fq.gz>: path to the truncated reverse gzipped FASTQ file\n" +
        "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)\n" +
        "\t\t<digest>: path to the digest file produced by the digest command\n" +
        "\t\t<active-digests>: path to a BED file with the coordinates of active digests\n" +
        String.format("\t\t<outfile>: optional name of output file (Default: \"%s.bam\")",DEFAULT_OUTPUT_BAM_NAME));
        System.out.println("\t\t<output-rejected>: output rejected reads to file (false if no -b option passed)\n");
    }

    private static void printHelpHeader() {
        String version=getVersion();
        System.out.println();
        System.out.println("Diachromatic (Analysis of Differential Capture Hi-C Interactions)\n"+
                "Version: "+version + "\n\n");
    }



    /**
     * Print usage information to provided OutputStream.
     */
    private static void printUsage(String message)
    {
        printHelpHeader();
        System.out.println(message + "\n\n"+
        "Usage: java -jar Diachromatic.jar <command> [options]\n\n"+
        "Available commands:\n\n");

        printDigestHelp(false);
        System.out.println();
        printTruncateHelp(false);
        System.out.println();
        printMapHelp(false);
        System.out.println();
        System.exit(0);
    }

    private void makeOutdirectoryIfNeeded() {
        File f = new File(outputDirectory);
        if (f.exists() && f.isFile()) {
            logger.error(String.format("Cannot make output directory called %s because a file of the same name exists", outputDirectory));
        } else if (! f.exists()) {
            f.mkdir(); // only make directory if necessary.
        }
    }
}
