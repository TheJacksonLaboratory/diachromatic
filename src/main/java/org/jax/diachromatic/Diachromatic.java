package org.jax.diachromatic;

import org.apache.log4j.Logger;
import org.jax.diachromatic.command.Command;
import org.jax.diachromatic.io.Commandline;

public class Diachromatic {
    static Logger logger = Logger.getLogger(Diachromatic.class.getName());


    public static void main(String args[]) {
        Diachromatic oncembobulator = new Diachromatic(args);
    }

    public Diachromatic(String args[]) {
        Commandline clp = new Commandline(args);
        Command command = clp.getCommand();
        command.execute();
    }
}
