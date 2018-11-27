package org.jax.diachromatic.command;

import org.jax.diachromatic.exception.DiachromaticException;

public abstract class Command {

    abstract public void execute() throws DiachromaticException;
}

