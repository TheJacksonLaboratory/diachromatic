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
            return (this.forwardDigest.getDigestEndPosition() == reverseDigest.getDigestStartPosition() - 1) || (this.reverseDigest.getDigestEndPosition() == forwardDigest.getDigestStartPosition() - 1);
        } else {
            return false;
        }
    }

    /**
     * @return true iff both digests in this digest pair are the same (artefact)
     */
    public boolean isSameDigest() {
        return this.forwardDigest.getChromosome().equals(this.reverseDigest.getChromosome()) &&
                this.forwardDigest.getDigestStartPosition() == this.reverseDigest.getDigestStartPosition() &&
                this.forwardDigest.getDigestEndPosition() == this.reverseDigest.getDigestEndPosition();
    }

    @Override
    public boolean equals(Object o) {
        if (o==null) return false;
        if (! (o instanceof DigestPair) ) return false;
        DigestPair other = (DigestPair) o;
        return (this.forwardDigest.equals(other.forwardDigest) && this.reverseDigest.equals(other.reverseDigest))
                ||
               (this.forwardDigest.equals(other.reverseDigest) && this.reverseDigest.equals(other.forwardDigest))
                ;
    }

    /** Hash code with lazily initialized value*/
    private int hashCode;
    @Override
    public int hashCode() {
        int result=hashCode;
        if (result==0) {
            if (forwardDigest.getDigestStartPosition() < reverseDigest.getDigestStartPosition()) {
                result = forwardDigest.hashCode();
                result = 31 * result + reverseDigest.hashCode();
                hashCode = result;
            } else {
                result = reverseDigest.hashCode();
                result = 31 * result + forwardDigest.hashCode();
                hashCode = result;
            }
        }
        return result;
    }

    public Digest get5digest() {
        if(forwardDigest.getDigestStartPosition()<reverseDigest.getDigestStartPosition()) {
            return forwardDigest;
        } else {
            return reverseDigest;
        }
    }

    public Digest get3digest() {
        if(forwardDigest.getDigestStartPosition()<reverseDigest.getDigestStartPosition()) {
            return reverseDigest;
        } else {
            return forwardDigest;
        }
    }


    @Override
    public String toString() {
        // fragment with the smaller starting position comes always first
        if(forwardDigest.getDigestStartPosition()<reverseDigest.getDigestStartPosition()) {
            return forwardDigest + "\t" + reverseDigest;
        } else {
            return reverseDigest + "\t" + forwardDigest;
        }
    }
}
