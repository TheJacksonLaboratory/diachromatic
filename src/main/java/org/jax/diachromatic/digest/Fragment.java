package org.jax.diachromatic.digest;

/**
 * This represents one fragment of a genomic digest.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @version 0.0.1
 */
public class Fragment implements Comparable<Fragment> {

    public final int enzymeNumber;
    public final int position;

    public Fragment(int enzymeNr, int pos) {
        enzymeNumber=enzymeNr;
        position=pos;
    }


    @Override
    public int compareTo(Fragment o) {
        return position - o.position;
    }
}
