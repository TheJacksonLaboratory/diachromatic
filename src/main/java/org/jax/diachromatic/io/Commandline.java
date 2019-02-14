package org.jax.diachromatic.io;

import java.io.File;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.AlignCommand;
import org.jax.diachromatic.command.CountCommand;
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

    /** Default size of margin of fragments used for calculating GC and repeat content. */
    private final static int DEFAULT_MARGIN_SIZE = 250; // also GC and repeat content will be calculated within GOPHER
    private final static String DEFAULT_OUTPUT_DIRECTORY="results";
    private final static String DEFAULT_FILENAME_PREFIX="prefix";

    /*
     options.addOption("od", "out-dir", true, "path to output directory") // general option
              // .addOption("op", "out-prefix", true, "outprefix for files in output directory") // general option
              // .addOption("h", "help", false, "shows help for current command") // general option

               .addOption("g", "genome", true, "path to genome FASTA file (with all chromosomes)") // digest specific option
               .addOption("m", "margin", true,"margin size for calculating GC and repeat content (default: 250 bp)") // digest specific option

               .addOption("e", "enzyme", true, "restriction enzyme name") // truncate specific option
                .addOption("s", "sticky-ends", false, "no fill-in of sticky ends was performed ") // truncate specific option
              / .addOption("q", "fastq-r1", true, "path to forward FASTQ input file") // truncate and align specific option
              / .addOption("r", "fastq-r2", true, "path to reverse FASTQ input file") // truncate and align specific option
               //.addOption("b", "bowtie-path", true, "path to bowtie2") // align specific option
               //.addOption("i", "bowtie-index", true, "path to bowtie2 index") // align specific option
             //  .addOption("bsu", "bowtie-stringent-unique", false, "use stringent settings for definition of uniquely mapped reads")
             //  .addOption("j", "bad", false, "output bad (rejected) reads to separated file") // align specific option
              // .addOption("d", "digest-file", true, "path to GOPHER digest file") // align (and count) specific option
               .addOption("o", "out", true, "name/path of output file/directory")
             //  .addOption("p", "thread-num", true, "number of threads used by bowtie2")
              // .addOption("l", "lower-frag-size-limit", true, "lower limit for fragment size")
              // .addOption("u", "upper-frag-size-limit", true, "upper limit for fragment size")
              // .addOption("v", "valid-pairs-bam", true, "path to BAM file with unique valid pairs produced using the 'align' subcommand")

     */

    @Deprecated //TODO--this was moved to GOPHER??
    private String enzyme=null;

//    private String pathToInputFastq1 = null;
//    private String pathToInputFastq2 = null;


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
//    private boolean outputRejectedReads=false;
//    /** path to the bowtie2 executable. */
//    private String bowtiepath=null;
//    private int marginsize;
//    private String pathToDiachromaticDigestFile=null;
//    private String pathToBowtieIndex=null;
//
//    private String pathToValidPairsBamFile=null;
//    private String pathToGopherDigestFile=null;
//
//    /**
//     * Lower and upper limit for size of valid pair fragments. Read pairs that correspond to hybrid fragments with
//     * a size outside this range will be categorized as wrong size pairs and discarded.
//     */
//    private Integer lowerFragSize = 50;
//    private Integer upperFragSize = 800;
//
//    private boolean useStringentUniqueSettings = false;

    public Commandline(String args[]) {
        /*
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
                printUsage("Error: Subcommand missing!");
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
            if (commandLine.hasOption("v")) {
                this.pathToValidPairsBamFile=commandLine.getOptionValue("v");
            }
            if (commandLine.hasOption("d")) {
                this.pathToDiachromaticDigestFile=commandLine.getOptionValue("d");
                this.pathToGopherDigestFile=commandLine.getOptionValue("d");
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
            if (commandLine.hasOption("bsu")) {
                this.useStringentUniqueSettings=true;
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
            outputDir = String.format("%s%s%s", outputDirectory, File.separator, filenamePrefix);
        }
        catch (ParseException parseException)  // checked exception
        {
            String msg = String.format("Could not parse options %s [%s]",clstring, parseException.toString());
           printUsage(msg );
        }
        if (doHelp) {
            switch (mycommand) {
                case "truncate": printHelpHeader(); printTruncateHelp(true); break;
                case "align": printHelpHeader(); printAlignHelp(true); break;
                case "count": printHelpHeader(); printCountHelp(true); break;
            }
            System.exit(1);
        }
        try {
            if (mycommand.equalsIgnoreCase("truncate")) {
                logger.trace(outputDir);
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
                this.command = new TruncateCommand(pathToInputFastq1, pathToInputFastq2, enzyme, stickyEnds, outputDir);

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
                        outputRejectedReads,
                        outputDir,
                        threadNum,
                        lowerFragSize,
                        upperFragSize,
                        filenamePrefix,
                        useStringentUniqueSettings
                        );
            } else if (mycommand.equalsIgnoreCase("count")) {
                if (this.pathToValidPairsBamFile == null) {
                    printCountHelp("ERROR: -v option required for count command. Please specify a BAM file with valid pairs.");
                } else if (this.pathToGopherDigestFile==null) {
                    printCountHelp("ERROR: -d option required for count command. Please specify a GOPHER digest file.");
                }
                logger.trace("List of arguments");
                logger.trace("=================");
                logger.trace("pathToValidPairsBamFile: " + this.pathToValidPairsBamFile);
                logger.trace("pathToGopherDigestFile: " + this.pathToGopherDigestFile);
                logger.trace("outputDir: " + this.outputDir);
                logger.trace("filenamePrefix: " + this.filenamePrefix + "\n");
                this.command=new CountCommand(this.pathToValidPairsBamFile,this.pathToGopherDigestFile,this.outputDir,this.filenamePrefix);
            } else {
                printUsage(String.format("Did not recognize command: %s", mycommand));
            }
        } catch (DiachromaticException de) {
            de.printStackTrace();
            printUsage(de.getMessage());
        }
         */
    }



    public Command getCommand() {
        return command;
    }

    /**
     * Construct and provide GNU-compatible Options.
     *
     * @return Options expected from command-line of GNU form.

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
               .addOption("bsu", "bowtie-stringent-unique", false, "use stringent settings for definition of uniquely mapped reads")
               .addOption("j", "bad", false, "output bad (rejected) reads to separated file") // align specific option
               .addOption("d", "digest-file", true, "path to GOPHER digest file") // align (and count) specific option
               .addOption("o", "out", true, "name/path of output file/directory")
               .addOption("p", "thread-num", true, "number of threads used by bowtie2")
               .addOption("l", "lower-frag-size-limit", true, "lower limit for fragment size")
               .addOption("u", "upper-frag-size-limit", true, "upper limit for fragment size")
               .addOption("v", "valid-pairs-bam", true, "path to BAM file with unique valid pairs produced using the 'align' subcommand")

        ;
        return options;
    } */

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
        System.out.println("Truncate:\n" +
        "\tjava -jar Diachromatic.jar truncate -q <forward.fq.gz> -r <reverse.fq.gz> -e <enzyme> \\ \n"+
            "\t\t\t[-sticky-ends] [-od <out-dir>] [-op <out-prefix>]\n\n"+

            "\t\t<forward.fq.gz>: path to the forward FASTQ file (may or may not be compressed with gzip)\n"+
            "\t\t<reverse.fq.gz>: path to the reverse FASTQ file (may or may not be compressed with gzip)\n"+
            "\t\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)\n"+
            "\t\t<true|false>: no fill-in of sticky ends was performed (Default: false)\n"+
            "\t\t<out-dir>: directory containing the output of the truncate command (Default: " + DEFAULT_OUTPUT_DIRECTORY + ")\n"+
            "\t\t<out-prefix>: prefix for all generated files in output directory (Default: " + DEFAULT_FILENAME_PREFIX + ")\n");
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
        System.out.println("Align:\n" +
        "\tjava -jar Diachromatic.jar align -b <bowtie2> -i <bowtie2-index> \\ \n" +

        "\t\t\t-q <forward.truncated.fq.gz> -r <reverse.truncated.fq.gz> \\ \n" +
        "\t\t\t-d <digest-file> [-od <outfile>] [-j <output-rejected>]\n" +
        "\t\t\t[-l <lower-frag-size-limit>] [-u <upper-frag-size-limit>]\n" +
        "\t\t\t[-o <outfile>] [-b] [-p] <thread-num>\n\n" +

        "\t\t<bowtie2>: path to bowtie2 executable\n" +
        "\t\t<bowtie2-index>: path to bowtie2 index of the reference genome\n" +
        "\t\t<bowtie-stringent-unique>: use stringent settings for definition of uniquely mapped reads\n" +
        "\t\t<forward.truncated.fq.gz>: path to the truncated forward gzipped FASTQ file\n" +
        "\t\t<reverse.truncated.fq.gz>: path to the truncated reverse gzipped FASTQ file\n" +
        "\t\t<lower-frag-size-limit>: lower limit for fragment size (Default: 150)\n" +
        "\t\t<upper-frag-size-limit>: upper limit for fragment size (Default: 800)\n" +
        "\t\t<digest-file>: path to the digest file produced with GOPHER\n" +
        "\t\t<thread-num>: number of threads used by bowtie2\n" +
        "\t\t<output-rejected>: output rejected reads to file)\n");
    }

    private static void printCountHelp(String message) {
        System.out.println("\n"+ message + "\n");
        printCountHelp(true);
        System.out.println();
        System.exit(0);
    }


    private static void printCountHelp(boolean specific) {
        if (specific) {
            System.out.println("The subcommand 'count' takes a BAM file containing unique valid read pairs determined during the align step and a digest file created with GOPHER and counts valid pairs between pairs of restriction fragments.\n");
        }
        System.out.println("Count:\n" +
                "\tjava -jar Diachromatic.jar count -v <valid-pairs-bam> -d <digest-file> \\ \n" +
                "\t\t\t[-od <out-dir>] [-op <out-prefix>]>\n\n" +

                "\t\t<valid-pairs-bam>: path to BAM file with unique valid read pairs produced during the align step (REQUIRED)\n" +
                "\t\t<digest-file>: path to the digest file produced using GOPHER (REQUIRED)\n" +
                "\t\t<out-dir>: directory for output files (Default: results)\n" +
                "\t\t<out-prefix>: prefix for names of output files (Default: prefix)\n");
    }

    private static void printHelpHeader() {
        String version=getVersion();
        System.out.println();
        System.out.println("Diachromatic (Analysis of Directed Capture Hi-C Interactions)\n"+
                "Version: "+version + "\n\n");
    }


    /**
     * Print usage information to provided OutputStream.
     */
    private static void printUsage(String message)
    {
        printHelpHeader();
        System.out.println(message + "\n\n"+
        "Usage: java -jar Diachromatic.jar <subcommand> [options]\n\n"+
        "Available commands:\n\n");

        printTruncateHelp(false);
        System.out.println();
        printAlignHelp(false);
        System.out.println();
        printCountHelp(false);
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
