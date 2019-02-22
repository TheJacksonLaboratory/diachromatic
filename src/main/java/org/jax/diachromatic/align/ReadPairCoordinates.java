package org.jax.diachromatic.align;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.Objects;

/**
 * This is a helper class of DeDupMap for the removal of duplicates that takes characteristics of Hi-C fragments
 * into account (orientation of read pairs matter).
 */
public class ReadPairCoordinates {
    private static final Logger logger = LogManager.getLogger();
    /**
     * Only the two 5' end positions and the orientation of read pairs are taken into account.
     */
    private final int fivePrimePos1;
    private final int fivePrimePos2;
    /** Code for orientation of read pairs. The code is set in the class {@link DeDupMap2} to a value between 0 and 4. */
    private final int readPairOrientation;
    /** Hash code (calculated in the constructor)*/
    private final int hashCode;


    ReadPairCoordinates(int fivePrimePos1, int fivePrimePos2, int readPairOrientation) {
        this.fivePrimePos1 = fivePrimePos1;
        this.fivePrimePos2 = fivePrimePos2;
        this.readPairOrientation = readPairOrientation;
        /* calculate the hashcode */
        int result=0;
        result=Integer.hashCode(readPairOrientation);
        result=31*result+Integer.hashCode(fivePrimePos1);
        result=31*result+Integer.hashCode(fivePrimePos2);
        this.hashCode=result;
    }

    /**
     * @param o the Object being compared with this.
     * @return true if o and this are equal
     */
    @Override
    public boolean equals(Object o) {
        if (o==null) return false;
        if (! (o instanceof ReadPairCoordinates) ) return false;
        ReadPairCoordinates other = (ReadPairCoordinates) o;
        if(this.readPairOrientation==other.readPairOrientation && this.fivePrimePos1==other.fivePrimePos1 && this.fivePrimePos2==other.fivePrimePos2)
        {
            /*
            logger.trace("false");
            logger.trace(this.readPairOrientation + "\t" + this.fivePrimePos1 + "\t" + this.fivePrimePos2);
            logger.trace(other.readPairOrientation + "\t" + other.fivePrimePos1 + "\t" + other.fivePrimePos2);
            */
            return true;
        } else {
            /*
            logger.trace("false");
            logger.trace(this.readPairOrientation + "\t" + this.fivePrimePos1 + "\t" + this.fivePrimePos2);
            logger.trace(other.readPairOrientation + "\t" + other.fivePrimePos1 + "\t" + other.fivePrimePos2);
            */
            return false;
        }
    }

    /**
     * {@link #hashCode} is calculate in the constructor {@link ReadPairCoordinates}.
     * @return hashcode for this object.
     */
    @Override
    public int hashCode() {
        return this.hashCode;
    }
}
