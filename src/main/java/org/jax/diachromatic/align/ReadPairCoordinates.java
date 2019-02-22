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
    Integer fivePrimePos1;
    Integer fivePrimePos2;
    Integer readPairOrientation;

    ReadPairCoordinates(int fivePrimePos1, int fivePrimePos2, int readPairOrientation) {
        this.fivePrimePos1 = fivePrimePos1;
        this.fivePrimePos2 = fivePrimePos2;
        this.readPairOrientation = readPairOrientation;
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
            return false;
        } else {
            /*
            logger.trace("false");
            logger.trace(this.readPairOrientation + "\t" + this.fivePrimePos1 + "\t" + this.fivePrimePos2);
            logger.trace(other.readPairOrientation + "\t" + other.fivePrimePos1 + "\t" + other.fivePrimePos2);
            */
            return true;
        }
    }

    /** Hash code with lazily initialized value*/
    private int hashCode;

    /**
     * Calculate the hash code based on the position of the Digest, the other values are not required
     * to obtain uniqueness.
     * @return hascode for this object.
     */
    @Override
    public int hashCode() {
        int result=hashCode;
        if (result==0) {
            result=readPairOrientation;
            result=31*result+fivePrimePos1;
            result=31*result+fivePrimePos2;
            hashCode=result;
            /*
            logger.trace("readPairOrientation: " + readPairOrientation);
            logger.trace("fivePrimePos1: " + fivePrimePos1);
            logger.trace("fivePrimePos2: " + fivePrimePos2);
            logger.trace("HashCode: " + result);
            */
        }
        return result;
    }
}
