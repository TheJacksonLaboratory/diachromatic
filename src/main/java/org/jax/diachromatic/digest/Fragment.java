package org.jax.diachromatic.digest;

/**
 * This represents one "cut" of a genomic digest, i.e., one position that is cut by a specific enzyme.
 * This can be used to generate a digest map which is an ordered list of cuts.
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
