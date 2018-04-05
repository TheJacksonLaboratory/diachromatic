package org.jax.diachromatic.map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import java.util.*;

/**
 * This class is intended for couting read pairs between interacting fragments.
 * for more than one condition.
 *
 * In essence this is a java HasMap with the fragment IDs as key and integer
 * arrays as values, but it has additional attributes that takes into account
 * more specific requirements.
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.1 (2018-04-04)
 */
public class InteractionCountsMap {
    private static final Logger logger = LogManager.getLogger();

    // fields
    // ------

    /**
     * Total number of conditions
     */
    private Integer number_of_conditions;

    /**
     * Hash map for counting interactions
     */
    private HashMap<String,List<Integer>> interaction_counts_map = null;

    /**
     * Number of interacting fragments
     */
    private Integer number_of_interacting_fragments;

    /**
     * Total current number of interactions for all conditions
     */
    private Integer interaction_count=0;


    // constructor
    // -----------

    public InteractionCountsMap(Integer number_of_conditions) {

        this.number_of_conditions = number_of_conditions;
        interaction_counts_map = new HashMap<String,List<Integer>>();
    }


    // private methods
    // ---------------

    /**
     * The interaction matrix is symmetric. To ensure that both fields of a given interaction are
     * counted together, the fragment with the smaller starting position comes always first.
     *
     * @param refID_1 Name of the reference sequence to which the first read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the first read is mapped.
     * @param refID_1 Name of the reference sequence to which the second read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the second read is mapped.
     *
     * @return Unique key for the given coordinates.
     */
    String getHashKey(String refID_1, Integer fragStaPos_1, String refID_2, Integer fragStaPos_2){

        String key="";
        String smallerRefID;
        String smallerFragPos;
        String largerRefID;
        String largerFragPos;

        // fragment with the smaller starting position comes always first
        if(fragStaPos_1 < fragStaPos_2) {
            smallerRefID=refID_1;largerRefID=refID_2;
            smallerFragPos=fragStaPos_1.toString();largerFragPos=fragStaPos_2.toString();
        } else {
            smallerRefID=refID_2;largerRefID=refID_1;
            smallerFragPos=fragStaPos_2.toString();largerFragPos=fragStaPos_1.toString();
        }

        // construct and return key
        key += smallerRefID;
        key += ":";
        key += smallerFragPos;
        key += ":";
        key += largerRefID;
        key += ":";
        key += largerFragPos;

        return key;
    }


    // public methods
    // --------------

    /**
     * This function is the main interface of this class.
     *
     * It takes as input:
     *
     * @param condition_num Identifier for the condition to be incremented.
     * @param refID_1 Name of the reference sequence to which the first read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the first read is mapped.
     * @param refID_1 Name of the reference sequence to which the second read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the second read is mapped.
     *
     * @throws Exception TODO: Add proper error handling for interactions within the same fragment. This should never happen (same internal artifact).
     *
     * @return The key for the incremented interaction is returned.
     *
     */
    public String incrementFragPair(Integer condition_num, String refID_1, Integer fragStaPos_1, String refID_2, Integer fragStaPos_2) {

        if(refID_1.compareTo(refID_2)==0 && fragStaPos_1==fragStaPos_2) {
            logger.warn("Interaction is within the same fragment. This should never happen for a valid pair (same internal artifact).");
            // TODO: Add proper error handling.
        }

        // generate unique key
        String hashKey = getHashKey(refID_1, fragStaPos_1, refID_2, fragStaPos_2);

        // check if hashKey exists
        if(!interaction_counts_map.containsKey(hashKey)) {

            // if not, init ArrayList of length that equals the total number of conditions with zero
            Integer[] integers = new Integer[number_of_conditions];
            Arrays.fill(integers, 0);
            List<Integer> newList = Arrays.asList(integers);
            interaction_counts_map.put(hashKey, newList);
        }

        // either way, increment associated array at corresponding position
        interaction_counts_map.get(hashKey).set(condition_num, interaction_counts_map.get(hashKey).get(condition_num)+1);
        interaction_count++;

        return hashKey;
    }

    /**
     *
     * @return Current total number of interactions for all conditions.
     */
    public Integer getCurrentTotalNumberOfInteractions() {
        return this.interaction_count;
    }

    public void printInteractionCountsMap() {

        Iterator it = interaction_counts_map.entrySet().iterator();

        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            System.out.println(pair.getKey() + " = " + pair.getValue());
            it.remove(); // avoids a ConcurrentModificationException
        }
    }

    /**
     * This function is primarily for testing and therefore it's public.
     *
     * @param hashKey
     * @param condition_id
     *
     * @return Current number of interactions for given fragment pair and condition.
     */
    public Integer getInteractionNumForKeyAndCondition(String hashKey, Integer condition_id) {
        Integer interNum = this.interaction_counts_map.get(hashKey).get(condition_id);
        return interNum;
        // TODO: Throws NullPointerException for some reasons. Check testGetInteractionNumForKeyAndCondition().
    }
}


