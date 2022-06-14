package org.jax.diachromatic.count;

public class BadReadCount {
    /** too short 5' */
    public int ts5count = 0;
    /** too short 3' */
    public int ts3count = 0;
    /** unligated 5' */
    public int ul5count = 0;
    /** unligated 3' */
    public int ul3count = 0;


    public void incrementTs5() {
        ts5count++;
    }

    public void incrementTs3() {
        ts3count++;
    }

    public void incrementUl5() {
        ul5count++;
    }

    public void incrementUl3() {
        ul3count++;
    }
}
