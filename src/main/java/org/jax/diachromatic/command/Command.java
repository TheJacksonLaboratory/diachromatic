package org.jax.diachromatic.command;

import com.beust.jcommander.Parameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.exception.DiachromaticException;

import java.io.File;

public abstract class Command {
    private static final Logger logger = LogManager.getLogger();

    @Parameter(names={"-x", "--prefix"},required = true,description = "prefix for files in output directory")
    protected String filenamePrefix;
    @Parameter(names={"-o", "--out-dir"},required = true,description = "name/path of output file/directory")
    protected String outputDir;

    abstract public void execute() throws DiachromaticException;

    public void makeOutdirectoryIfNeeded() {
        File f = new File(outputDir);
        if (f.exists() && f.isFile()) {
            logger.error(String.format("Cannot make output directory called %s because a file of the same name exists", outputDir));
        } else if (! f.exists()) {
            f.mkdir(); // only make directory if necessary.
        }
    }
}

