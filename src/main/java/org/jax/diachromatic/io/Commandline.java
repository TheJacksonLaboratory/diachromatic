package org.jax.diachromatic.io;

import java.io.PrintWriter;
import java.io.Writer;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.apache.commons.cli.*;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.command.DigestCommand;
import org.jax.diachromatic.command.MapCommand;
import org.jax.diachromatic.command.TruncateCommand;
import org.jax.diachromatic.exception.DiachromaticException;

/**
 * Class to capture options and command from the command line.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.2 (2018-01-05)
 */
public class Commandline {
    private static final Logger logger = LogManager.getLogger();
    private Command command = null;
    /**
     * The default name of the file that is produced by the {@code digest} command.
     */
    private final static String DEFAULT_DIGEST_FILE_NAME = "diachromaticDigest.txt";

    private final static String DEFAULT_OUTPUT_DIRECTORY = "results";

    private final static String DEFAULT_TRUNCATION_SUFFIX = "truncated";

    private final static String DEFAULT_OUTPUT_BAM_NAME = "diachromatic-processed";


    private String genomeDirectory = null;
    private String outputFilePath = null;
    private String outputDirectory = null;
    private String enzyme = null;
    /**
     * path to the bowtie2 executable.
     */
    private String bowtiepath = null;
    private String pathToBowtieIndex = null;
    private String pathToInputFastq1 = null;
    private String pathToInputFastq2 = null;
    private String pathToDiachromaticDigestFile = null;
    private String suffix = null;
    private boolean outputRejectedReads = false;

    public Commandline(String args[]) {
        final CommandLineParser cmdLineGnuParser = new DefaultParser();

        final Options gnuOptions = constructGnuOptions();
        org.apache.commons.cli.CommandLine commandLine;

        String mycommand = null;
        String clstring = "";
        if (args != null && args.length > 0) {
            clstring = Arrays.stream(args).collect(Collectors.joining(" "));
        }
        try {
            commandLine = cmdLineGnuParser.parse(gnuOptions, args);
            String category[] = commandLine.getArgs();
            if (category.length == 0) {
                printUsage("command missing");
            } else if (category.length > 1) {
                String cmd = Arrays.stream(category).collect(Collectors.joining("\t"));
                System.out.println(String.format("%d arguments for command: %s", category.length, cmd));
            } else {
                mycommand = category[0];

            }
            if (commandLine.getArgs().length < 1) {
                printUsage("no arguments passed");
                return;
            }
            if (commandLine.hasOption("g")) {
                this.genomeDirectory = commandLine.getOptionValue("g");
            }
            if (commandLine.hasOption("o")) {
                this.outputFilePath = commandLine.getOptionValue("o");
            }
            if (commandLine.hasOption("outdir")) {
                outputDirectory = commandLine.getOptionValue("outdir");
            }
            if (commandLine.hasOption("e")) {
                this.enzyme = commandLine.getOptionValue("e");
            }
            if (commandLine.hasOption("b")) {
                this.bowtiepath = commandLine.getOptionValue("b");
            }
            if (commandLine.hasOption("d")) {
                this.pathToDiachromaticDigestFile = commandLine.getOptionValue("d");
            }
            if (commandLine.hasOption("b")) {
                outputRejectedReads = true;
            } else {
                outputRejectedReads = false;
            }
            if (commandLine.hasOption("i")) {
                this.pathToBowtieIndex = commandLine.getOptionValue("i");
            }
            if (commandLine.hasOption("q")) {
                this.pathToInputFastq1 = commandLine.getOptionValue("q");
            }
            if (commandLine.hasOption("r")) {
                this.pathToInputFastq2 = commandLine.getOptionValue("r");
            }
            if (commandLine.hasOption("s")) {
                this.suffix = commandLine.getOptionValue("s");
            }
        } catch (ParseException parseException)  // checked exception
        {
            String msg = String.format("Could not parse options %s [%s]", clstring, parseException.toString());
            printUsage(msg);
        }
        try {
            if (mycommand.equals("digest")) {
                if (this.genomeDirectory == null) {
                    printUsageSubprogram("-g option required for digest command", "digest");
                }
                if (this.enzyme == null) {
                    printUsageSubprogram("-e option required for digest command", "digest");
                }
                if (this.outputFilePath == null) {
                    outputFilePath = DEFAULT_DIGEST_FILE_NAME;
                }
                this.command = new DigestCommand(this.genomeDirectory, enzyme, this.outputFilePath);

            } else if (mycommand.equalsIgnoreCase("truncate")) {
                if (this.outputDirectory == null) {
                    this.outputDirectory = DEFAULT_OUTPUT_DIRECTORY;
                } else if (this.pathToInputFastq1 == null) {
                    printUsageSubprogram("-q option required for truncate command", "truncate");
                } else if (this.pathToInputFastq2 == null) {
                    printUsageSubprogram("-r option required for truncate command", "truncate");
                } else if (this.enzyme == null) {
                    printUsageSubprogram("-e option required for truncate command", "truncate");
                }
                if (suffix == null) {
                    suffix = DEFAULT_TRUNCATION_SUFFIX;
                }
                //String outdir, String file1, String file2, String enzymeName
                this.command = new TruncateCommand(outputDirectory, pathToInputFastq1, pathToInputFastq2, enzyme, suffix);
            } else if (mycommand.equalsIgnoreCase("map")) {
                logger.trace("IN MAP");
                if (this.bowtiepath == null) {
                    printUsageSubprogram("-b option required for map command", "map");
                }
                if (this.pathToBowtieIndex == null) {
                    printUsageSubprogram("-i option (bowtie index) required for map command", "map");
                }
                if (this.pathToInputFastq1 == null) {
                    printUsageSubprogram("-q option (FASTQ 1) required for map command", "map");
                }
                if (this.pathToInputFastq2 == null) {
                    printUsageSubprogram("-r option (FASTQ 2) required for map command", "map");
                }
                if (this.outputFilePath == null) {
                    outputFilePath = DEFAULT_OUTPUT_BAM_NAME;
                }
                if (pathToDiachromaticDigestFile == null) {
                    printUsageSubprogram("-d option required for map command", "map");
                }
                this.command = new MapCommand(bowtiepath,
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
    public static Options constructGnuOptions() {
        final Options gnuOptions = new Options();
        gnuOptions.addOption("o", "out", true, "name/path of output file/directory")
                .addOption("g", "genome", true, "genome directory (with FASTA files)")
                .addOption("e", "enzyme", true, "restriction enzyme name")
                .addOption("b", "bowtie", true, "path to bowtie2")
                .addOption("d", "digest", true, "path to diachromatic digest file")
                .addOption("s", "suffix", true, "suffix for output filenames")
                .addOption("i", "bowtieindex", true, "path to bowtie2 index")
                .addOption("outdir", "outdir", true, "path to output directory")
                .addOption("q", "q", true, "path to forward FASTQ input file")
                .addOption("r", "r", true, "path to reverse FASTQ input file")
                .addOption("j", "bad", false, "output bad (reJected) reads to separated file")
                .addOption(Option.builder("f1").longOpt("file1").desc("path to fastq file 1").hasArg(true).argName("file1").build())
                .addOption(Option.builder("f2").longOpt("file2").desc("path to fastq file 2").hasArg(true).argName("file2").build());
        return gnuOptions;
    }

    public static String getVersion() {
        String version = "0.0.0";// default, should be overwritten by the following.
        try {
            Package p = Commandline.class.getPackage();
            version = p.getImplementationVersion();
        } catch (Exception e) {
            // do nothing
        }
        return version;
    }


    private static void printMap(PrintWriter writer) {
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
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s.bam\")", DEFAULT_OUTPUT_BAM_NAME));
        writer.println("\t-b: output rejected reads to file (false if no -b option passed)");
        writer.println();
    }

    private static void printTruncate(PrintWriter writer) {
        writer.println("truncate:");
        writer.println("\tjava -jar Diachromatic.jar truncate -q <forward.fq.gz> \\");
        writer.println("\t\t-r <reverse.fq.gz> -e <enzyme> -s <suffix> --outdir <directory>");
        writer.println("\t<forward.fq.gz>: path to the forward FASTQ file (may or may not be compressed with gzip)");
        writer.println("\t<reverse.fq.gz>: path to the reverse FASTQ file (may or may not be compressed with gzip)");
        writer.println("\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)");
        writer.println("\t<suffix>: suffix that will be added to the output truncated FASTQ files");
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s\")", DEFAULT_TRUNCATION_SUFFIX));
        writer.println();
    }

    private static void printDigest(PrintWriter writer) {
        writer.println("digest:");
        writer.println("\tjava -jar Diachromatic.jar digest -g <path> -e <enzyme> [-o <outfile>]");
        writer.println("\t<path>: path to a directory containing indexed genome FASTA files");
        writer.println("\t<enzyme>: symbol of the restriction enzyme (e.g., DpnII)");
        writer.println(String.format("\t<outfile>: optional name of output file (Default: \"%s\")", DEFAULT_DIGEST_FILE_NAME));
        writer.println();
    }


    /**
     * Print usage information to provided OutputStream.
     */
    public static void printUsage(String message) {
        String version = getVersion();
        final PrintWriter writer = new PrintWriter(System.out);
        writer.println(message);
        writer.println();
        //usageFormatter.printUsage(writer, 120, applicationName, options);
        writer.println("Program: Diachromatic (Analysis of Differential Capture Hi-C Interactions)");
        writer.println("Version: " + version);
        writer.println();
        writer.println("Usage: java -jar Diachromatic.jar <command> [options]");
        writer.println();
        writer.println("Available commands:");
        writer.println();
        printDigest(writer);
        printTruncate(writer);
        printMap(writer);
        writer.close();
        System.exit(0);
    }

    public static void printUsageSubprogram(String message, String subprogramname) {
        String version = getVersion();
        final PrintWriter writer = new PrintWriter(System.out);
        writer.println();
        //usageFormatter.printUsage(writer, 120, applicationName, options);
        writer.println("Program: Diachromatic (Analysis of Differential Capture Hi-C Interactions)");
        writer.println("Version: " + version);
        writer.println();
        writer.println("[ERROR]\n\t" + message);
        writer.println();
        switch (subprogramname) {
            case "truncate":
                printTruncate(writer);
            case "digest":
                printDigest(writer);
            case "map":
                printMap(writer);
            default:
                writer.println("Did not recognize command: " + subprogramname);
        }
        writer.close();
        System.exit(0);
    }

}
