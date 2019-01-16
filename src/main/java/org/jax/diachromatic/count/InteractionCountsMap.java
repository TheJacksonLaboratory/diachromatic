package org.jax.diachromatic.count;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jax.diachromatic.exception.IncrementSameInternalInteractionException;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * This class is intended for counting read pairs between interacting fragments for one or more conditions.
 *
 * In essence this is a java HasMap with keys assembled from the coordinates of the interacting fragments
 * and Integer arrays that contain the numbers of interactions for individual conditions,
 * but it has additional features that take into account more specific requirements.
 *
 * Once the counting of interactions is done, a second hash align can be optionally derived for the numbers
 * of reads at interacting fragments. The class provides methods that can be used to write the
 * content of the respective hash maps to a text file. Furthermore, numbers such as the total number
 * of interactions for a given condition or the total number of interacting fragments for a given condition
 * are determined and can be queried via public methods.
 *
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
     * Total number of read pair orientations
     */
    private Integer number_of_orientations = 4;

    /**
     * Counter reads within active fragments
     */
    private Integer[] read_count = null;

    /**
     * Counter reads within active fragments
     */
    private Integer[] active_read_count = null;

    /**
     * Total current number of interactions for each condition
     */
    private Integer[] interaction_count = null;

    /**
     * Total current number of singleton interactions for each condition.
     * Note: Initialized only after execution of the function
     * 'printInteractionCountsMapAsCountTable', which is not nice.
     */
    private Integer[] n_singleton_interactions = null;

    /**
     * Number of interactions between two active fragments for each condition
     */
    private Integer[] active_active_interaction_count = null;

    /**
     * Number of interactions between two inactive fragments for each condition
     */
    private Integer[] inactive_inactive_interaction_count = null;

    /**
     * Number of interactions between active and inactive fragments for each condition (both directions)
     */
    private Integer[] active_inactive_interaction_count = null;

    /**
     * Total number of interacting fragments for each condition
     */
    private Integer[] interacting_fragment_count = null;

    /**
     * Total number of active interacting fragments for each condition
     */
    private Integer[] active_interacting_fragment_count = null;

    /**
     * Hash align for counting interactions
     */
    private HashMap<String,List<Integer>> interaction_counts_map = null;

    /**
     * Space efficient hash for counting interactions using pairs of Integer keys. For a given key pair,
     * the smaller key points to a HashMap and within this HashMap the larger key points to an count array.
     * Interactions are counted separately for different read pair orientations.
     */
    private HashMap<Integer,HashMap<Integer,short[]>> interaction_key_counts_map = null;

    /**
     * Output filename for interaction counts
     */
    //private String interactionCountsTableFileName = "diachromatic.interaction.counts.table.tsv";



    /**
     * Hash align for read counts at interacting fragments
     */
    private HashMap<String,List<Integer>> fragment_interaction_counts_map = new HashMap<String,List<Integer>>();

    /**
     * Output filename read counts at interacting fragments
     */
    //private String interactingFragmentsCountsTableFileName = "diachromatic.interacting.fragments.counts.table.tsv";



    // constructor
    // -----------

    public InteractionCountsMap(Integer number_of_conditions) {

        this.number_of_conditions = number_of_conditions;
        interaction_counts_map = new HashMap<String,List<Integer>>();

        interaction_key_counts_map = new HashMap<Integer,HashMap<Integer,short[]>>();

        this.interaction_count = new Integer[number_of_conditions];
        Arrays.fill(interaction_count, 0);

        this.n_singleton_interactions = new Integer[number_of_conditions];
        Arrays.fill(n_singleton_interactions, 0);

        this.active_active_interaction_count = new Integer[number_of_conditions];
        Arrays.fill(active_active_interaction_count, 0);

        this.inactive_inactive_interaction_count = new Integer[number_of_conditions];
        Arrays.fill(inactive_inactive_interaction_count, 0);

        this.active_inactive_interaction_count = new Integer[number_of_conditions];
        Arrays.fill(active_inactive_interaction_count, 0);

        this.interacting_fragment_count = new Integer[number_of_conditions];
        Arrays.fill(interacting_fragment_count, 0);

        this.active_interacting_fragment_count = new Integer[number_of_conditions];
        Arrays.fill(active_interacting_fragment_count, 0);

        this.active_read_count = new Integer[number_of_conditions];
        Arrays.fill(active_read_count, 0);

        this.read_count = new Integer[number_of_conditions];
        Arrays.fill(read_count, 0);
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
     * @return Unique key for the given coordinates of a fragment pair.
     */
    String getHashKey(String refID_1, Integer fragStaPos_1, Integer fragEndPos_1, boolean fragActive_1, String refID_2, Integer fragStaPos_2, Integer fragEndPos_2, boolean fragActive_2){

        String key="";
        String smallerRefID;
        String smallerFragStaPos;
        String smallerFragEndPos;
        String smallerActivationState;
        String largerRefID;
        String largerFragStaPos;
        String largerFragEndPos;
        String largerActivationState;

        // fragment with the smaller starting position comes always first
        if(fragStaPos_1 < fragStaPos_2) {
            smallerRefID=refID_1;
            smallerFragStaPos=fragStaPos_1.toString();
            smallerFragEndPos=fragEndPos_1.toString();
            if(fragActive_1) {smallerActivationState="A";} else {smallerActivationState="I";}
            largerRefID=refID_2;
            largerFragStaPos=fragStaPos_2.toString();
            largerFragEndPos=fragEndPos_2.toString();
            if(fragActive_2) {largerActivationState="A";} else {largerActivationState="I";}

        } else {
            smallerRefID=refID_2;
            smallerFragStaPos=fragStaPos_2.toString();
            smallerFragEndPos=fragEndPos_2.toString();
            if(fragActive_2) {smallerActivationState="A";} else {smallerActivationState="I";}
            largerRefID=refID_1;
            largerFragStaPos=fragStaPos_1.toString();
            largerFragEndPos=fragEndPos_1.toString();
            if(fragActive_1) {largerActivationState="A";} else {largerActivationState="I";}
        }

        // construct and return key
        key += smallerRefID; key += ":"; key += smallerFragStaPos; key += "-"; key += smallerFragEndPos; key += ":" ; key += smallerActivationState; key += ";"; key += largerRefID; key += ":"; key += largerFragStaPos; key += "-"; key += largerFragEndPos;key += ":" ; key += largerActivationState;
        return key;
    }

    private boolean fragmentKeyIsActive(String fragmentKey) {
        String tmp[] = fragmentKey.split(":");
        if(tmp[2].equals("A")) {
            return true;
        } else {
            return false;
        }
    }


    // public methods
    // --------------

    public void incrementFragPair2(Integer condition_num, Integer digestKey1, Integer digestKey2, int relOriTag) {

        // ensure that (a,b) and (b,a) pairs are treated as the same
        Integer firstKey, secondKey;
        if (digestKey1 < digestKey2) {
            firstKey = digestKey1;
            secondKey = digestKey2;
        } else {
            firstKey = digestKey2;
            secondKey = digestKey1;
        }

        // create new structures, if necessary
        HashMap<Integer, short[]> newInnerMap;
        short[] newCountArray;
        if (!interaction_key_counts_map.containsKey(firstKey)) {
            // create new inner hash map
            newInnerMap = new HashMap<Integer, short[]>();
            // create new count array and put on inner hash map
            newCountArray = new short[4];
            //Arrays.fill(newCountArray, 0);
            newInnerMap.put(secondKey, newCountArray);
            // and put inner hash on outer hash
            interaction_key_counts_map.put(firstKey, newInnerMap);
            interaction_count[condition_num]++;
        } else if (!interaction_key_counts_map.get(firstKey).containsKey(secondKey)) {
            // there is already an inner hash map only count array needs to be created
            newCountArray = new short[4];
            //Arrays.fill(newCountArray, 0);
            // put count array on inner hash map
            interaction_key_counts_map.get(firstKey).put(secondKey,newCountArray);
            interaction_count[condition_num]++;
        }

        // increment either way
        interaction_key_counts_map.get(firstKey).get(secondKey)[relOriTag]++;

        //logger.trace(digestKey1 + "," + digestKey2 + " -> " + interaction_key_counts_map.get(firstKey).get(secondKey)[relOriTag]);
    }





    /**
     * This method is the main interface of this class. It can be used to increment the count of interactions
     * for a given pair of restriction fragment and condition.
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
     * TODO: Try to pass the corresponding two digests to this function instead of the long list of arguments.
     *
     */
    public String incrementFragPair(Integer condition_num, String refID_1, Integer fragStaPos_1, Integer fragEndPos_1, boolean fragActive_1, String refID_2, Integer fragStaPos_2, Integer fragEndPos_2, boolean fragActive_2, String relOriTag) throws IncrementSameInternalInteractionException {

        // generate unique String key
        String hashKey = getHashKey(refID_1, fragStaPos_1, fragEndPos_1, fragActive_1, refID_2, fragStaPos_2, fragEndPos_2, fragActive_2);

        // count interaction separately for simple and twisted
        String oriTag;
        if(relOriTag.equals("F1F2") || relOriTag.equals("F2F1") || relOriTag.equals("R1R2") || relOriTag.equals("R2R1")) {
            oriTag="T"; // twisted loop
        } else {
            oriTag="S"; // simple loop
        }

        hashKey += ";";
        hashKey += oriTag;

        // count reads in active fragments
        if(fragActive_1) {active_read_count[condition_num]++;}
        if(fragActive_2) {active_read_count[condition_num]++;}
        read_count[condition_num]=read_count[condition_num]+2;


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

            // if this the first read pair observed for the current condition, also increment total interaction count for current condition
            if(interaction_counts_map.get(hashKey).get(condition_num).intValue()==1) {
                interaction_count[condition_num]++;
                String frag[] = hashKey.split(";");
                if (fragmentKeyIsActive(frag[0]) && fragmentKeyIsActive(frag[1])) {
                    active_active_interaction_count[condition_num]++;
                } else if (!fragmentKeyIsActive(frag[0]) && !fragmentKeyIsActive(frag[1])) {
                    inactive_inactive_interaction_count[condition_num]++;
                } else {
                    active_inactive_interaction_count[condition_num]++;
                }
            }
        }
        catch (IncrementSameInternalInteractionException e) {
            logger.warn("IncrementSameInternalInteraction occurred. Interaction is within the same fragment.");
        }
        return hashKey;
    }

    /**
     *
     * @return Total number of interactions for a given condition.
     *
     */
    public Integer getTotalNumberOfInteractionsForCondition(Integer condition) {
        return this.interaction_count[condition];
    }

    /**
     *
     * @return Total number of interactions with only one read pair for a given condition.
     *
     */
    public Integer getTotalNumberOfSingletonInteractionsForCondition(Integer condition) {
        return this.n_singleton_interactions[condition];
    }

    /**
     *
     * @return Total number of interactions between two active fragments for a given condition.
     *
     */
    public Integer getNumberOfInteractionsBetweenActiveFragmentsForCondition(Integer condition) {
        return this.active_active_interaction_count[condition];
    }

    /**
     *
     * @return Total number of interactions between two active fragments for a given condition.
     *
     */
    public Integer getNumberOfInteractionsBetweenInactiveFragmentsForCondition(Integer condition) {
        return this.inactive_inactive_interaction_count[condition];
    }

    /**
     *
     * @return Total number of interactions between active and inactive fragments for a given condition (both direction).
     *
     */
    public Integer getNumberOfInteractionsBetweenActiveAndInactiveFragmentsForCondition(Integer condition) {
        return this.active_inactive_interaction_count[condition];
    }

    /**
     *
     * @return Total number of interacting fragments for a given condition.
     *
     */
    public Integer getTotalNumberOfInteractingFragmentsForCondition(Integer condition) {
        return this.interacting_fragment_count[condition];
    }

    /**
     *
     * @return Total number of active interacting fragments for a given condition.
     *
     */
    public Integer getTotalNumberOfActiveInteractingFragmentsForCondition(Integer condition) {
        return this.active_interacting_fragment_count[condition];
    }

    public double getTargetEnrichmentCoefficientForCondition(Integer condition) {
        return active_read_count[condition]/read_count[condition].doubleValue();
    }

    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a pair of interacting fragments.
     * The first six fields contain the coordinates of the two interacting fragments.
     * All following fields contain the number of read pairs observed for the given fragment pair
     * for individual conditions.
     *
     */
    public void printInteractionCountsMapAsCountTable(String interactionCountsTableFileName) throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(interactionCountsTableFileName));

        // iterate over hash and write to file
        Iterator it = interaction_counts_map.entrySet().iterator();
        HashSet seenHashKeys = new HashSet(); // keep track of seen interaction, either simple or twisted, to avoid double printing
        while (it.hasNext()) {

            Map.Entry pair = (Map.Entry)it.next();

            String hashKey = pair.getKey().toString();
            String[] frags = hashKey.split(";");
            String baseHashKey = frags[0];
            baseHashKey += ";";
            baseHashKey += frags[1];

            if(seenHashKeys.contains(baseHashKey)) {
                continue;
            } else {
                //seenHashKeys.add(baseHashKey);
            }

            String[] tmp1 = frags[0].split(":");
            String refID_1 = tmp1[0];
            String[] tmp2 = tmp1[1].split("-");
            String fragStaPos_1 = tmp2[0];
            String fragEndPos_1 = tmp2[1];
            String fragActivationState_1 = tmp1[2];

            tmp1 = frags[1].split(":");
            String refID_2 = tmp1[0];
            tmp2 = tmp1[1].split("-");
            String fragStaPos_2 = tmp2[0];
            String fragEndPos_2 = tmp2[1];
            String fragActivationState_2 = tmp1[2];

            String relOriTag = frags[2];

            printStream.print(refID_1 + "\t" + fragStaPos_1 + "\t" + fragEndPos_1 + "\t" + fragActivationState_1 + "\t");
            printStream.print(refID_2 + "\t" + fragStaPos_2 + "\t" + fragEndPos_2 + "\t" + fragActivationState_2);

            for(int i=0; i<number_of_conditions; i++) {
                Integer simpleNum;
                Integer twistedNum;
                if(relOriTag.equals("S")) {
                    simpleNum = interaction_counts_map.get(hashKey).get(i);
                    // check if there is also a twisted interaction
                    String twistedHaskKey = baseHashKey.concat(";T");
                    if(interaction_counts_map.containsKey(twistedHaskKey)) {
                        twistedNum = interaction_counts_map.get(twistedHaskKey).get(i);
                    } else {
                        twistedNum = 0;
                    }
                } else {
                    twistedNum = interaction_counts_map.get(hashKey).get(i);
                    // check if there is also a twisted interaction
                    String simpleHaskKey = baseHashKey.concat(";S");
                    if(interaction_counts_map.containsKey(simpleHaskKey)) {
                        simpleNum = interaction_counts_map.get(simpleHaskKey).get(i);
                    } else {
                        simpleNum = 0;
                    }
                }

                printStream.print("\t");
                printStream.print(simpleNum);
                printStream.print(":");
                printStream.print(twistedNum);

                if(simpleNum+twistedNum==1) {n_singleton_interactions[i]++;}

                //printStream.print(interaction_counts_map.get(hashKey).get(i));
            }
            printStream.print("\n");

        }
    }

    public void printInteractionCountsMapAsCountTable2(DigestMap digestMap, String interactionCountsTableFileName) throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(interactionCountsTableFileName));

        // iterate over all interactions
        for (Integer key1 : interaction_key_counts_map.keySet()) {
            for(Integer key2 : interaction_key_counts_map.get(key1).keySet()) {
                //logger.trace(key1 + ", " + key2);
                //logger.trace(digestMap.getCoordinatesForDigestKey(key1) + ", " + digestMap.getCoordinatesForDigestKey(key2));
                //logger.trace("Interaction count: " + interaction_key_counts_map.get(key1).get(key2)[0]);
                printStream.print(digestMap.getCoordinatesForDigestKey(key1) + "\t" + digestMap.getCoordinatesForDigestKey(key2) + "\t");
                short [] countArray = interaction_key_counts_map.get(key1).get(key2);
                int simple=0, twisted=0;
                for(int i=0; i<countArray.length; i++) {
                    printStream.print(countArray[i] + "\t");
                    if(i<2) {
                        twisted=twisted+countArray[i];
                    } else {
                        simple=simple+countArray[i];
                    }
                }
                printStream.print(simple + ":" + twisted + "\n");

            }
        }
    }

    /**
     * This method derives read counts at interacting fragments from the interaction counts.
     */
    public void deriveReadCountsAtInteractingFragments() {

        int cnt=0;
        Iterator it = interaction_counts_map.entrySet().iterator();
        while (it.hasNext()) {
            cnt++;
            if(cnt%1000==0) {logger.trace(cnt);}
            Map.Entry pair = (Map.Entry)it.next();
            String hashKey = pair.getKey().toString();
            String[] fragKey = hashKey.split(";");

            // check if fragment keys exists and create array if necessary
            if(!fragment_interaction_counts_map.containsKey(fragKey[0])) {
                Integer[] frag1counts = new Integer[number_of_conditions];
                Arrays.fill(frag1counts, 0);
                List<Integer> newList = Arrays.asList(frag1counts);
                fragment_interaction_counts_map.put(fragKey[0], newList);
            }
            if(!fragment_interaction_counts_map.containsKey(fragKey[1])) {
                Integer[] frag1counts = new Integer[number_of_conditions];
                Arrays.fill(frag1counts, 0);
                List<Integer> newList = Arrays.asList(frag1counts);
                fragment_interaction_counts_map.put(fragKey[1], newList);
            }

            // increment counts for both fragments and each condition
            for(int i = 0; i<number_of_conditions; i++) {
                if(fragment_interaction_counts_map.get(fragKey[0]).get(i)==0) {
                    interacting_fragment_count[i]++;
                    if(fragmentKeyIsActive(fragKey[0])) {
                        active_interacting_fragment_count[i]++;
                    }
                }
                Integer newVal = fragment_interaction_counts_map.get(fragKey[0]).get(i) + interaction_counts_map.get(hashKey).get(i);
                fragment_interaction_counts_map.get(fragKey[0]).set(i,newVal);
                if(fragment_interaction_counts_map.get(fragKey[1]).get(i)==0) {
                    interacting_fragment_count[i]++;
                    if(fragmentKeyIsActive(fragKey[1])) {
                        active_interacting_fragment_count[i]++;
                    }
                }
                newVal = fragment_interaction_counts_map.get(fragKey[1]).get(i) + interaction_counts_map.get(hashKey).get(i);
                fragment_interaction_counts_map.get(fragKey[1]).set(i,newVal);
            }
        }
    }


    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a restriction fragment that contains at least
     * one read of a valid pair for at least one condition.
     *
     * The first three fields contain the coordinates of the interacting fragment.
     * All following fields contain the number of reads (that are part of a valid pair)
     * at the interacting fragments individual conditions.
     */
    public void printFragmentInteractionCountsMapAsCountTable(String interactingFragmentsCountsTableFileName) throws FileNotFoundException {

        // derive counts
        this.deriveReadCountsAtInteractingFragments();

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(interactingFragmentsCountsTableFileName));

        Iterator it = fragment_interaction_counts_map.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            String hashKey = pair.getKey().toString();
            String[] tmp = hashKey.split(":");
            String[] tmp2 = tmp[1].split("-");
             printStream.print(tmp[0] + "\t" + tmp2[0] + "\t" + tmp2[1] + "\t" + tmp[2]);
            for(int j = 0; j<number_of_conditions; j++) {
                printStream.print("\t" + fragment_interaction_counts_map.get(hashKey).get(j));
            }
            printStream.print("\n");
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
    public Integer getNumberOfInteractionsForKeyAndCondition(String hashKey, Integer condition_id) {
        Integer interNum = interaction_counts_map.get(hashKey).get(condition_id);
        return interNum;
    }

    /**
     * This method is primarily for testing and therefore it's public.
     *
     * @param hashKey
     * @param condition_id
     *
     * @return Current number of reads for a given fragment and condition.
     */
    public Integer getNumberOfReadsAtInteractingFragmentForKeyAndCondition(String hashKey, Integer condition_id) {
        Integer interNum = fragment_interaction_counts_map.get(hashKey).get(condition_id);
        return interNum;
    }
}
