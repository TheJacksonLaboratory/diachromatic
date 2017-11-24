package org.jax.diachromatic;


import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.io.Commandline;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
/**
 * An application to process Hi-C data for differential reads counts in fragments surrounding the
 * transcription start site using probe design by VPV.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.1 (2017-11-15)
 */
public class Diachromatic {
    private static final Logger logger = LogManager.getLogger();

    public static void main(String args[]) {
        Commandline clp = new Commandline(args);
        Command command = clp.getCommand();
        logger.trace(String.format("running command %s",command));
        command.execute();
    }
}
