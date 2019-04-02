package org.jax.diachromatic.util;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A simple status bar that only work on terminals where "\r" has an affect.
 *
 * The progress is done/shown in the closed interval <code>[min, max]</code>.
 *
 * @author <a href="mailto:manuel.holtgrewe@charite.de">Manuel Holtgrewe</a>
 */
public final class ProgressBar {

    private static final Logger logger = LogManager.getLogger();
    /** smallest value */
    private final long min;
    /** largest value */
    private final long max;

    /** Initialize progress bar with the given settings */
    public ProgressBar(long min, long max) {
        this.min = min;
        this.max = max;

    }

    /** @return smallest value to represent */
    public long getMin() {
        return min;
    }

    /** @return largest value to represent */
    public long getMax() {
        return max;
    }


    /** print progress up to position <code>pos</code>. */
    public void print(long pos) {
        int percent = (int) Math.ceil(100.0 * (pos - this.min) / (this.max - this.min));
        StringBuilder bar = new StringBuilder("[");

        for (int i = 0; i < 50; i++) {
            if (i < (percent / 2)) {
                bar.append("=");
            } else if (i == (percent / 2)) {
                bar.append(">");
            } else {
                bar.append(" ");
            }
        }

        bar.append("]   " + percent + "%     ");
        logger.error("\r" + bar.toString());
        if (pos == max)
            logger.error("\n");
    }

    public void finish() {
        logger.error("\n");
    }

}
