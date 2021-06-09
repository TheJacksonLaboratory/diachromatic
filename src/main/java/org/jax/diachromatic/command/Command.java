package org.jax.diachromatic.command;


import java.io.File;

import picocli.CommandLine;

public class Command {

    @CommandLine.Option(names = {"-o", "--out-dir"}, required = true, description = "Name/path of summarize directory. Will be created, if not yet existing.")
    protected String outputDir = "";

    @CommandLine.Option(names = {"-x", "--prefix"}, required = true, description = "Prefix for files in summarize directory.")
    protected String filenamePrefix;


    public void makeOutdirectoryIfNeeded() {
        File f = new File(outputDir);
        if (f.exists() && f.isFile()) {
            System.err.printf("Cannot make summarize directory called %s because a file of the same name exists", outputDir);
        } else if (!f.exists()) {
            f.mkdir(); // only make directory if necessary.
        }
    }
}
