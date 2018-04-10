package org.jax.diachromatic.map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jax.diachromatic.exception.IncrementSameInternalInteractionException;

import javax.rmi.CORBA.Util;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * This class is intended for counting read pairs between interacting fragments
 * for one or more conditions.
 *
 * In essence this is a java HasMap with keys assembled from the coordinates of the interacting fragments
 * and Integer arrays that contain the numbers of interactions, but it has additional features that take into account
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


    private String interactionCountsTableFileName = "diachromatic.interaction.counts.table.tsv";


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
     * This method assembles the key accordingly.
     *
     * @param refID_1 Name of the reference sequence to which the first read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the first read is mapped.
     * @param fragEndPos_1 Last position of the fragment to which the first read is mapped.
     * @param refID_2 Name of the reference sequence to which the second read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_2 Starting position of the fragment to which the second read is mapped.
     * @param fragEndPos_2 Last position of the fragment to which the first read is mapped.
     *
     * @return Unique key for the given coordinates.
     */
    String getHashKey(String refID_1, Integer fragStaPos_1, Integer fragEndPos_1, String refID_2, Integer fragStaPos_2, Integer fragEndPos_2){

        String key="";
        String smallerRefID;
        String smallerFragStaPos;
        String smallerFragEndPos;
        String largerRefID;
        String largerFragStaPos;
        String largerFragEndPos;

        // fragment with the smaller starting position comes always first
        if(fragStaPos_1 < fragStaPos_2) {
            smallerRefID=refID_1;
            smallerFragStaPos=fragStaPos_1.toString();
            smallerFragEndPos=fragEndPos_1.toString();
            largerRefID=refID_2;
            largerFragStaPos=fragStaPos_2.toString();
            largerFragEndPos=fragEndPos_2.toString();
        } else {
            smallerRefID=refID_2;
            smallerFragStaPos=fragStaPos_2.toString();
            smallerFragEndPos=fragEndPos_2.toString();
            largerRefID=refID_1;
            largerFragStaPos=fragStaPos_1.toString();
            largerFragEndPos=fragEndPos_1.toString();
        }

        // construct and return key
        key += smallerRefID; key += ":"; key += smallerFragStaPos; key += "-"; key += smallerFragEndPos; key += ";"; key += largerRefID; key += ":"; key += largerFragStaPos; key += "-"; key += largerFragEndPos;
        return key;
    }


    // public methods
    // --------------

    /**
     * This method is the main interface of this class.
     *
     * It takes as input:
     *
     * @param condition_num Identifier for the condition to be incremented.
     * @param refID_1 Name of the reference sequence to which the first read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the first read is mapped.
     * @param refID_1 Name of the reference sequence to which the second read of the pair is mapped, usually the name of a chromosome, e.g. chr1.
     * @param fragStaPos_1 Starting position of the fragment to which the second read is mapped.
     *
     * @throws IncrementSameInternalInteractionException
     *
     * @return The key for the incremented interaction is returned.
     *
     */
    public String incrementFragPair(Integer condition_num, String refID_1, Integer fragStaPos_1, Integer fragEndPos_1, String refID_2, Integer fragStaPos_2, Integer fragEndPos_2) throws IncrementSameInternalInteractionException {

        // generate unique key
        String hashKey = getHashKey(refID_1, fragStaPos_1, fragEndPos_1, refID_2, fragStaPos_2, fragEndPos_2);

        try {

            if(refID_1.compareTo(refID_2)==0 && fragStaPos_1==fragStaPos_2) {
                throw new IncrementSameInternalInteractionException();
            }

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
        }
        catch (IncrementSameInternalInteractionException e) {
            logger.warn("IncrementSameInternalInteraction occured. Interaction is within the same fragment.");
        }
        return hashKey;
    }

    /**
     *
     * @return Current total number of interactions for all conditions.
     *
     */
    public Integer getCurrentTotalNumberOfInteractions() {
        return this.interaction_count;
    }

    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a pair of interacting fragments.
     * The first six fields contain the coordinates of the two interacting fragments.
     * All following fields contain the number of read pairs observed for the given fragment pair
     * for individual conditions.
     *
     * TODO: Include information about active/inactive fragments as soon as this information is available.
     *
     */
    public void printInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(interactionCountsTableFileName));

        // iterate over hash and write to file
        Iterator it = interaction_counts_map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            String hashKey = pair.getKey().toString();
            String[] frags = hashKey.split(";");

            String[] tmp1 = frags[0].split(":");
            String refID_1 = tmp1[0];
            String[] tmp2 = tmp1[1].split("-");
            String fragStaPos_1 = tmp2[0];
            String fragEndPos_1 = tmp2[1];

            tmp1 = frags[1].split(":");
            String refID_2 = tmp1[0];
            tmp2 = tmp1[1].split("-");
            String fragStaPos_2 = tmp2[0];
            String fragEndPos_2 = tmp2[1];

            printStream.print(refID_1 + "\t" + fragStaPos_1 + "\t" + fragEndPos_1 + "\t");
            printStream.print(refID_2 + "\t" + fragStaPos_2 + "\t" + fragEndPos_2);

            for(int i=0; i<number_of_conditions; i++) {
                printStream.print("\t");
                printStream.print(interaction_counts_map.get(hashKey).get(i));
            }
            printStream.print("\n");
        }
    }

    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a restriction fragment that contains at least
     * one read of a valid pair for at least one condition.
     * The first three fields contain the coordinates of the interacting fragment.
     * All following fields contain the number of reads (that are part of a valid pair)
     * that are mapped to the fragment for individual conditions.
     *
     */
    public void printFragmentInteractionCountsMapAsCountTable() throws FileNotFoundException {

        // sort hash keys lexicographically
        Collection<String> keySet = interaction_counts_map.keySet();
        List sortedKeyList = new ArrayList(keySet);
        Collections.sort(sortedKeyList);

        String[] tmp = String.valueOf(sortedKeyList.get(0)).split(";");
        String prev_frag=tmp[0];
        System.out.println(prev_frag);
        for(int i=1; i<sortedKeyList.size(); i++) {
            String pairKey = String.valueOf(sortedKeyList.get(i));
            tmp = pairKey.split(";");
            String curr_frag = tmp[0];
            System.out.print(curr_frag + "\t");
            for(int j = 0; j<number_of_conditions;j++) {
                System.out.print("\t" + interaction_counts_map.get(pairKey).get(j));
            }
            System.out.print("\n");
        }




    }



    /**
     * This method is primarily for testing and therefore it's public.
     *
     * @param hashKey
     * @param condition_id
     *
     * @return Current number of interactions for given fragment pair and condition.
     */
    public Integer getInteractionNumForKeyAndCondition(String hashKey, Integer condition_id) {
        Integer interNum = interaction_counts_map.get(hashKey).get(condition_id);
        return interNum;
    }

}
