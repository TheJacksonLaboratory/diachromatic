package org.jax.diachromatic.map;

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

    // fields
    // ------

    /**
     * Number of conditions
     */
    private Integer number_of_conditions;

    /**
     * Hash for counting interactions
     */
    private HashMap<String,List<Integer>> interaction_counts_map = null;

    /**
     * Number of interacting fragments
     */
    private Integer number_of_interacting_fragments;


    // constructor
    // -----------

    public InteractionCountsMap(Integer number_of_conditions) {

        this.number_of_conditions = number_of_conditions;
        interaction_counts_map = new HashMap<String,List<Integer>>();
    }


    // private methods
    // ---------------

    String getHashKey(String refID_1, Integer fragNum_1, String refID_2, Integer fragNum_2){

        String key="";
        String smallerRefID;
        String smallerFragNum;
        String largerRefID;
        String largerFragNum;

        // lexicographically smaller reference ID first
        if(refID_1.compareTo(refID_2) <= 0) {
            smallerRefID=refID_1;largerRefID=refID_2;
        } else {
            smallerRefID=refID_2;largerRefID=refID_1;
        }

        // numerically smaller fragment number first
        if(fragNum_1 < fragNum_2) {
            smallerFragNum=fragNum_1.toString();largerFragNum=fragNum_2.toString();
        } else {
            smallerFragNum=fragNum_2.toString();largerFragNum=fragNum_1.toString();
        }

        // construct and return key
        key += smallerRefID;
        key += ":";
        key += smallerFragNum;
        key += ":";
        key += largerRefID;
        key += ":";
        key += largerFragNum;
        return key;
    }


    // public methods
    // --------------

    /**
     * This is the main interface of this class.
     *
     * It returns the key for the incremented interaction.
     *
     */
    public String incrementFragPair(Integer condition_num, String refID_1, Integer fragNum_1, String refID_2, Integer fragNum_2) {

        // generate unique key
        String hashKey = getHashKey(refID_1, fragNum_1, refID_2, fragNum_2);

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
        return hashKey;
    }

    public void printInteractionCountsMap() {

        Iterator it = interaction_counts_map.entrySet().iterator();

        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            System.out.println(pair.getKey() + " = " + pair.getValue());
            it.remove(); // avoids a ConcurrentModificationException
        }
    }
}
