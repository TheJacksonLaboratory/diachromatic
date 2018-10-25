package org.jax.diachromatic.io;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.apache.commons.cli.*;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.DigestCommand;
import org.jax.diachromatic.command.AlignCommand;
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
    private final static String DEFAULT_DIGEST_FILE_NAME = "default"; // digestion and fasta indexing will be moved to GOPHER
    /** Default size of margin of fragments used for calculating GC and repeat content. */
    private final static int DEFAULT_MARGIN_SIZE = 250; // also GC and repeat content will be calculated within GOPHER
    /** Absolute path to the combined genome fasta file (which will be indexed only if necessary). */
    private String genomeFastaFile=null; // only needed for digestion

    private final static String DEFAULT_OUTPUT_DIRECTORY="results";
    private final static String DEFAULT_FILENAME_PREFIX="prefix";


    private String enzyme=null;

    private String pathToInputFastq1 = null;
    private String pathToInputFastq2 = null;


    /**
     * general fields
     */
    private boolean doHelp=false;
    private String outputFilePath = null;

    private String filenamePrefix = DEFAULT_FILENAME_PREFIX;
    private String outputDirectory = DEFAULT_OUTPUT_DIRECTORY;
    private String outputPathPrefix = null;
    private Integer threadNum=1;
    private boolean stickyEnds=false;

    private boolean filledEnds = true; //

    /**
     * align specific fields
     */
    private boolean outputRejectedReads=false;
    /** path to the bowtie2 executable. */
    private String bowtiepath=null;
    private int marginsize;
    private String pathToDiachromaticDigestFile=null;
    private String pathToActiveDigestsFile=null;
    private String pathToBowtieIndex=null;

    /**
     * Lower and upper limit for size of valid pair fragments. Read pairs that correspond to hybrid fragments with
     * a size outside this range will be categorized as wrong size pairs and discarded.
     */
    private Integer lowerFragSize = 150;
    private Integer upperFragSize = 800;

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
            if (commandLine.hasOption("od")) {
                this.outputDirectory=commandLine.getOptionValue("od");
            }
            if (commandLine.hasOption("op")) {
                this.filenamePrefix=commandLine.getOptionValue("op");
            }

            if (commandLine.hasOption("b")) {
                this.bowtiepath=commandLine.getOptionValue("b");
            }
            if (commandLine.hasOption("p")) {
                this.threadNum=Integer.parseInt(commandLine.getOptionValue("p"));
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
            if (commandLine.hasOption("l")) {
                this.lowerFragSize=Integer.parseInt(commandLine.getOptionValue("l"));
            }
            if (commandLine.hasOption("u")) {
                this.upperFragSize=Integer.parseInt(commandLine.getOptionValue("u"));
            }
            if (commandLine.hasOption("s")) {
                this.stickyEnds=true;
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
                    printUsage("[ERROR] -m option requires integer value. You passed \"" + m + "\"");
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
            if (commandLine.hasOption("filenamePrefix")) {
                this.filenamePrefix=commandLine.getOptionValue("filenamePrefix");
            }

            // create directory for output in any case
            makeOutdirectoryIfNeeded();

            // create prefix for output files including the path to the output directory
            outputPathPrefix = String.format("%s%s%s", outputDirectory, File.separator, filenamePrefix);
        }
        catch (ParseException parseException)  // checked exception
        {
            String msg = String.format("Could not parse options %s [%s]",clstring, parseException.toString());
           printUsage(msg );
        }
        if (doHelp) {
            switch (mycommand) {
                case "digest": printHelpHeader(); printDigestHelp(true); break;
                case "truncate": printHelpHeader(); printTruncateHelp(true); break;
                case "align": printHelpHeader(); printAlignHelp(true); break;
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
                logger.trace(outputPathPrefix);
                if (this.outputDirectory == null) {
                    this.outputDirectory=DEFAULT_OUTPUT_DIRECTORY;
                } else if (this.pathToInputFastq1 == null) {
                    printTruncateHelp("-q option required for truncate command");
                } else if (this.pathToInputFastq2 == null) {
                    printTruncateHelp("-r option required for truncate command");
                } else if (this.enzyme == null) {
                    printTruncateHelp("-e option required for truncate command");
                }
                if (filenamePrefix==null) {
                    filenamePrefix=DEFAULT_FILENAME_PREFIX;
                }
                this.command = new TruncateCommand(pathToInputFastq1, pathToInputFastq2, enzyme, stickyEnds, outputPathPrefix);

            } else if (mycommand.equalsIgnoreCase("align")) {
                if (this.bowtiepath == null) {
                    printAlignHelp("-b option required for align command");
                }
                if (this.pathToBowtieIndex == null) {
                    printAlignHelp("-i option (bowtie index) required for align command");
                }
                if (this.pathToInputFastq1 == null) {
                    printAlignHelp("-q option (FASTQ 1) required for align command");
                }
                if (this.pathToInputFastq2 == null) {
                    printAlignHelp("-r option (FASTQ 2) required for align command");
                }
                if (pathToDiachromaticDigestFile == null) {
                    printAlignHelp("-d option required for align command");
                }
                if(threadNum<0 || threadNum>100) {
                    printAlignHelp("number of threads has to be an integer between 1 and 100");
                }
                this.command=new AlignCommand(
                        bowtiepath,
                        pathToBowtieIndex,
                        pathToInputFastq1,
                        pathToInputFastq2,
                        pathToDiachromaticDigestFile,
                        pathToActiveDigestsFile,
                        outputRejectedReads,
                        outputPathPrefix,
                        threadNum,
                        lowerFragSize,
                        upperFragSize,
                        filenamePrefix
                        );
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
        options.addOption("od", "out-dir", true, "path to output directory") // general option
               .addOption("op", "out-prefix", true, "outprefix for files in output directory") // general option
               .addOption("h", "help", false, "shows help for current command") // general option

               .addOption("g", "genome", true, "path to genome FASTA file (with all chromosomes)") // digest specific option
               .addOption("m", "margin", true,"margin size for calculating GC and repeat content (default: 250 bp)") // digest specific option

               .addOption("e", "enzyme", true, "restriction enzyme name") // truncate specific option
                .addOption("s", "sticky-ends", false, "no fill-in of sticky ends was performed ") // truncate specific option
               .addOption("q", "fastq-r1", true, "path to forward FASTQ input file") // truncate and align specific option
               .addOption("r", "fastq-r2", true, "path to reverse FASTQ input file") // truncate and align specific option
               .addOption("b", "bowtie-path", true, "path to bowtie2") // align specific option
               .addOption("i", "bowtie-index", true, "path to bowtie2 index") // align specific option
               .addOption("j", "bad", false, "output bad (rejected) reads to separated file") // align specific option
               .addOption("d", "digest", true, "path to GOPHER digest file") // align (and count) specific option
               .addOption("a", "active-digests", true, "path to BED file with active digests (overwrites information in digest file)") // align (and count) specific option
               .addOption("o", "out", true, "name/path of output file/directory")
               .addOption("p", "thread-num", true, "number of threads used by bowtie2")
               .addOption("l", "lower-frag-size-limit", true, "lower limit for fragment size")
               .addOption("u", "upper-frag-size-limit", true, "upper limit for fragment size")

        ;
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
            "\t\t\t-r <reverse.fq.gz> -e <enzyme> [-out-dir <output_directory>] [-out-prefix <filename_prefix>]\n\n"+

            "\t\t<forward.fq.gz>: path to the forward FASTQ file (may or may not be compressed with gzip)\n"+
            "\t\t<reverse.fq.gz>: path to the reverse FASTQ file (may or may not be compressed with gzip)\n"+
            "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)\n"+
            "\t\t<output_directory>: directory containing the output of the truncate command (Default: " + DEFAULT_OUTPUT_DIRECTORY + ")\n"+
            "\t\t<filename_prefix>: prefix for all generated files in output directory (Default: " + DEFAULT_FILENAME_PREFIX + ")\n");
    }

    private static void printAlignHelp(String message) {
        System.out.println("\n"+ message + "\n");
        printAlignHelp(true);
        System.out.println();
        System.exit(0);
    }


    private static void printAlignHelp(boolean specific) {
        if (specific) {
            System.out.println("Command and options for align function");
            System.out.println("\talign align uses bowtie2 to align reads and then performs Q/C and repairing.");
        }
        System.out.println("\n" +
        "\tjava -jar Diachromatic.jar align -b <bowtie2> -i <bowtie2-index> \\ \n" +

        "\t\t\t-q <forward.truncated.fq.gz> -r <reverse.truncated.fq.gz> \\ \n" +
        "\t\t\t-d <digest> [-od <outfile>] [-j <output-rejected>]\n" +
        "\t\t\t[-l <lower-frag-size-limit>] [-u <upper-frag-size-limit>]\n" +
        "\t\t\t[-a <active-digests>] [-o <outfile>] [-b] [-p] <thread-num>\n\n" +

        "\t\t<bowtie2>: path to bowtie2 executable\n" +
        "\t\t<bowtie2-index>: path to bowtie2 index for digested genome\n" +
        "\t\t<forward.truncated.fq.gz>: path to the truncated forward gzipped FASTQ file\n" +
        "\t\t<reverse.truncated.fq.gz>: path to the truncated reverse gzipped FASTQ file\n" +
        "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII or HindIII)\n" +
        "\t\t<lower-frag-size-limit>: lower limit for fragment size (Default: 150)\n" +
        "\t\t<upper-frag-size-limit>: upper limit for fragment size (Default: 800)\n" +
        "\t\t<digest>: path to the digest file produced by the digest command\n" +
        "\t\t<active-digests>: path to a BED file with the coordinates of active digests\n" +
        "\t\t<thread-num>: number of threads used by bowtie2\n" +
        "\t\t<output-rejected>: output rejected reads to file)\n");
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
        printAlignHelp(false);
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
