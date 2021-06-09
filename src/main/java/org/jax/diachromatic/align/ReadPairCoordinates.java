package org.jax.diachromatic.align;

/**
 * This is a helper class of DeDupMap for the removal of duplicates that takes characteristics of Hi-C fragments
 * into account (orientation of read pairs matter).
 */
public class ReadPairCoordinates {
    /**
     * Only the two 5' end positions and the orientation of read pairs are taken into account.
     */
    private final int fivePrimePos1;
    private final int fivePrimePos2;
    /** Code for orientation of read pairs. The code is set in the class {@link DeDupMap} to a value between 0 and 4. */
    private final int readPairOrientation;
    /** Hash code (calculated in the constructor)*/
    private final int hashCode;


    ReadPairCoordinates(int fivePrimePos1, int fivePrimePos2, int readPairOrientation) {
        this.fivePrimePos1 = fivePrimePos1;
        this.fivePrimePos2 = fivePrimePos2;
        this.readPairOrientation = readPairOrientation;
        /* calculate the hashcode */
        int result;
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
        return (this.readPairOrientation==other.readPairOrientation
                && this.fivePrimePos1==other.fivePrimePos1
                && this.fivePrimePos2==other.fivePrimePos2);
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
