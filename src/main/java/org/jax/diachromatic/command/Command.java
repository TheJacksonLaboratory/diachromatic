package org.jax.diachromatic.command;

import com.beust.jcommander.Parameter;
import org.jax.diachromatic.exception.DiachromaticException;

public abstract class Command {
    @Parameter(names={"-x", "--prefix"},required = true,description = "prefix for files in output directory")
    protected String filenamePrefix;
    @Parameter(names={"-o", "--out-dir"},required = true,description = "name/path of output file/directory")
    protected String outputPath;

    abstract public void execute() throws DiachromaticException;
}

