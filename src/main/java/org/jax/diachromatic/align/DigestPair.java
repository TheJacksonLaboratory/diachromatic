package org.jax.diachromatic.align;


/**
 * In capture Hi-C, each valid read pair is mapped to two restriction digests. This class contains the Digests
 * that align with the forward and reverse reads of a readpair.
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.1.2 (2018-01-06)
 */
public class DigestPair {
    /** The digest that matches with the forward read. */
    private final Digest forwardDigest;
    /** The digest that matches with the reverse read. */
    private final Digest reverseDigest;


    public DigestPair(Digest forward, Digest reverse) {
        forwardDigest=forward;
        reverseDigest=reverse;
    }

    public Digest forward() { return forwardDigest;}

    public Digest reverse() { return reverseDigest; }



    // int max_possible_insert_size = digestPair.forward().getSize() + digestPair.reverse().getSize();
    public int getMaximumPossibleInsertSize() {
        return forwardDigest.getSize() + reverseDigest.getSize();
    }

    public boolean isAdjacent() {
        if(this.forwardDigest.getChromosome().equals(this.reverseDigest.getChromosome())) {
            if((this.forwardDigest.getEndpos() == reverseDigest.getStartpos()-1)) {// || (this.reverseDigest.getEndpos() == forwardDigest.getStartpos()-1)) {
                return true;
            } else {
                return  false;
            }
        } else {
            return false;
        }
    }
}
