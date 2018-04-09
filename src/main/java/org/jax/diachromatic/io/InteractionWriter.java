package org.jax.diachromatic.io;

import org.jax.diachromatic.map.InteractionCountsMap;

import java.util.Iterator;
import java.util.Map;

public class InteractionWriter {

    /**
     * This function takes an object of class InteractionCountsMap and writes for each condition
     * one text file in iBED format to disk.
     * See: http://regulatorygenomicsgroup.org/wp-content/uploads/Chicago_vignette.html#output-files
     * For a description of the iBED format.
     *
     *
     */
    public void writeInteractionsInIbedFormat(InteractionCountsMap inputMap) {

    }

    /**
     * This function takes an object of class InteractionCountsMap and writes for each condition
     * one text file in seqmonk format to disk.
     * See: http://regulatorygenomicsgroup.org/wp-content/uploads/Chicago_vignette.html#output-files
     * For a description of the iBED format.
     *
     *
     */
    public void writeInteractionsInSeqMonkFormat(InteractionCountsMap inputMap) {

    }

    /**
     * This function takes an object of class InteractionCountsMap and writes one tab delimited file
     * text format to disk.
     *
     * The first six columns contain the coordinates of the two interacting fragments.
     *
     *
     */
    public void writeInteractionsAsInteractionTable(InteractionCountsMap inputMap) {


    }

}
