package org.jax.diachromatic.command;

import org.jax.diachromatic.Diachromatic;
import org.jax.diachromatic.exception.DiachromaticException;
import org.jax.diachromatic.map.Bowtie2Runner;

public class MapCommand extends Command {

    private final String bowtiepath;

    public MapCommand(String bowtie){
        bowtiepath=bowtie;

    }

    public void execute() {
        try {
            Bowtie2Runner runner = new Bowtie2Runner(bowtiepath);
            runner.run();
        } catch (DiachromaticException e){
            e.printStackTrace();
        }
    }
    @Override
    public String toString() {return "diachromatic:map";}
}
