package org.jax.diachromatic.command;


import org.jax.diachromatic.digest.RestrictionEnzyme;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.truncation.Truncator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import picocli.CommandLine;

import java.io.File;
import java.util.List;
import java.util.concurrent.Callable;

import static org.jax.diachromatic.digest.RestrictionEnzyme.parseRestrictionEnzymes;

/**
 * Class to coordinate truncation of chimeric Hi-C reads containing ligation junctions.

 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.2 (2018-01-05)
 */

@CommandLine.Command(name = "truncate",
        aliases = {"T"},
        mixinStandardHelpOptions = true,
        description = "Truncate reads at the positions of Hi-C religation sequences.")
public class TruncateCommand extends Command implements Callable<Integer> {
    private static final Logger logger = LoggerFactory.getLogger(TruncateCommand.class);

    @CommandLine.Option(names={"-q","fastq-r1"}, required = true, description = "Path to forward FASTQ input file.", order = 3)
    private String fastaqFile1;
    @CommandLine.Option(names={"-r","fastq-r2"}, required = true, description = "Path to reverse FASTQ input file.", order = 4)
    private String fastaqFile2;
    @CommandLine.Option(names={"-e", "--enzyme"}, required = true, description = "Restriction enzyme name.", order = 5)
    private String enzymeName;
    @CommandLine.Option(names={"-s", "--sticky-ends"},description = "No fill-in of sticky ends was performed.", order = 6)
    private boolean stickyEnds=false;

    private Truncator truncator = null;
    private RestrictionEnzyme re = null;

    public TruncateCommand(){}

    private void init() throws DiachromaticException {
        List<RestrictionEnzyme>  enzymelist = parseRestrictionEnzymes();
        re=enzymelist.stream().filter(r->r.getName().equalsIgnoreCase(enzymeName)).findFirst().orElse(null);
        if (re==null) {
            throw new DiachromaticException(String.format("Could not identify restriction enzyme for \"%s\"",enzymeName));
        }
        File f = new File(fastaqFile1);
        if(!f.exists()) {
            throw new DiachromaticException(String.format("%s does not exist", fastaqFile1));
        }
        f = new File(fastaqFile2);
        if(!f.exists()) {
            throw new DiachromaticException(String.format("%s does not exist", fastaqFile2));
        }
        String outputDirAndFilePrefix=String.format("%s%s%s", outputDir, File.separator,filenamePrefix);
        truncator = new Truncator(fastaqFile1,fastaqFile2, re, stickyEnds, outputDirAndFilePrefix);
    }


    @Override
    public Integer call() {
        makeOutdirectoryIfNeeded();
        logger.trace(String.format("Starting truncate command on files %s and %s",fastaqFile1,fastaqFile2));
        logger.trace(outputDir);
        try {
            init();
            truncator.parseFASTQ();
        } catch (DiachromaticException e) {
            logger.error("Error encountered while truncating FASTQ reads: {}", e.getLocalizedMessage());
        }
        return 0;
    }

    @Override
    public String toString() {return "diachromatic:truncate";}
}
