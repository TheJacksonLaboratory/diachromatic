package org.jax.diachromatic.count;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.jax.diachromatic.exception.IncrementSameInternalInteractionException;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * This class is intended for counting read pairs between interacting fragments.
 *
 * In essence this is a java HasMap with keys assembled from the coordinates of the interacting fragments
 * and Integer arrays that contain the numbers of interactions,
 * but it has additional features that take into account more specific requirements.
 *
 * Once the counting of interactions is done, a second hash align can be optionally derived for the numbers
 * of reads at interacting fragments. The class provides methods that can be used to write the
 * content of the respective hash maps to a text file. Furthermore, numbers such as the total number
 * of interactions or the total number of interacting fragments
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
     * Counter reads within active fragments
     */
    private Integer read_count = null;

    /**
     * Counter reads within active fragments
     */
    private Integer active_read_count = null;

    /**
     * Total current number of interactions
     */
    private Integer interaction_count = null;

    /**
     * Total current number of singleton interactions
     * Note: Initialized only after execution of the function
     * 'printInteractionCountsMapAsCountTable'.
     */
    private Integer n_singleton_interactions = null;

    /**
     * Number of interactions between two active fragments
     */
    private Integer active_active_interaction_count = null;

    /**
     * Number of interactions between two inactive fragments
     */
    private Integer inactive_inactive_interaction_count = null;

    /**
     * Number of interactions between active and inactive fragments (both directions)
     */
    private Integer active_inactive_interaction_count = null;

    /**
     * Total number of interacting fragments
     */
    private Integer interacting_fragment_count = null;

    /**
     * Total number of active interacting fragments
     */
    private Integer active_interacting_fragment_count = null;

    /**
     * Hash align for counting interactions
     */
    private HashMap<String,Integer> interaction_counts_map = null;

    /**
     * Output filename for interaction counts
     */
    //private String interactionCountsTableFileName = "diachromatic.interaction.counts.table.tsv";



    /**
     * Hash align for read counts at interacting fragments
     */
    private HashMap<String,Integer> fragment_interaction_counts_map = new HashMap<String,Integer>();

    /**
     * Output filename read counts at interacting fragments
     */
    //private String interactingFragmentsCountsTableFileName = "diachromatic.interacting.fragments.counts.table.tsv";



    // constructor
    // -----------

    public InteractionCountsMap() {

        interaction_counts_map = new HashMap<String,Integer>();
        this.interaction_count = 0;
        this.n_singleton_interactions = 0;
        this.active_active_interaction_count = 0;
        this.inactive_inactive_interaction_count = 0;
        this.active_inactive_interaction_count = 0;
        this.interacting_fragment_count = 0;
        this.active_interacting_fragment_count = 0;
        this.active_read_count = 0;
        this.read_count = 0;
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
        String[] tmp = fragmentKey.split(":");
        return tmp[2].equals("A");
    }


    // public methods
    // --------------

    /**
     * This method is the main interface of this class. It can be used to increment the count of interactions
     * for a given pair of restriction fragments.
     *
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
    public String incrementFragPair(String refID_1, Integer fragStaPos_1, Integer fragEndPos_1, boolean fragActive_1, String refID_2, Integer fragStaPos_2, Integer fragEndPos_2, boolean fragActive_2, String relOriTag) throws IncrementSameInternalInteractionException {

        // generate unique key
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
        if(fragActive_1) {active_read_count++;}
        if(fragActive_2) {active_read_count++;}
        read_count = read_count+2;


        try {

            if(refID_1.compareTo(refID_2)==0 && fragStaPos_1==fragStaPos_2) {
                throw new IncrementSameInternalInteractionException();
            }

           // check if hashKey exists
            if(!interaction_counts_map.containsKey(hashKey)) {

                // if not, init with zero
                interaction_counts_map.put(hashKey, 0);
            }

            // either way, increment associated array at corresponding position
            interaction_counts_map.put(hashKey,interaction_counts_map.get(hashKey)+1);


            // if this the first read pair observed, also increment total interaction count
            if(interaction_counts_map.get(hashKey)==1) {
                interaction_count++;
                String[] frag = hashKey.split(";");
                if (fragmentKeyIsActive(frag[0]) && fragmentKeyIsActive(frag[1])) {
                    active_active_interaction_count++;
                } else if (!fragmentKeyIsActive(frag[0]) && !fragmentKeyIsActive(frag[1])) {
                    inactive_inactive_interaction_count++;
                } else {
                    active_inactive_interaction_count++;
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
     * @return Total number of interactions
     *
     */
    public Integer getTotalNumberOfInteractions() {
        return this.interaction_count;
    }

    /**
     *
     * @return Total number of interactions with only one read pair.
     *
     */
    public Integer getTotalNumberOfSingletonInteractions() {
        return this.n_singleton_interactions;
    }

    /**
     *
     * @return Total number of interactions between two active fragments.
     *
     */
    public Integer getNumberOfInteractionsBetweenActiveFragments() {
        return this.active_active_interaction_count;
    }

    /**
     *
     * @return Total number of interactions between two active fragments.
     *
     */
    public Integer getNumberOfInteractionsBetweenInactiveFragments() {
        return this.inactive_inactive_interaction_count;
    }

    /**
     *
     * @return Total number of interactions between active and inactive fragments (both direction).
     *
     */
    public Integer getNumberOfInteractionsBetweenActiveAndInactiveFragments() {
        return this.active_inactive_interaction_count;
    }

    /**
     *
     * @return Total number of interacting fragments.
     *
     */
    public Integer getTotalNumberOfInteractingFragments() {
        return this.interacting_fragment_count;
    }

    /**
     *
     * @return Total number of active interacting fragments.
     *
     */
    public Integer getTotalNumberOfActiveInteractingFragments() {
        return this.active_interacting_fragment_count;
    }

    public double getTargetEnrichmentCoefficient() {
        return active_read_count/read_count.doubleValue();
    }

    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a pair of interacting fragments.
     * The first six fields contain the coordinates of the two interacting fragments.
     * All following fields contain the number of read pairs observed for the given fragment pair.
     *
     */
    public void printInteractionCountsMapAsCountTable(String interactionCountsTableFileName) throws FileNotFoundException {

        // create file for output
        PrintStream printStream = new PrintStream(new FileOutputStream(interactionCountsTableFileName));

        // iterate over hash and write to file
        Iterator it = interaction_counts_map.entrySet().iterator();
        HashSet<String> seenHashKeys = new HashSet(); // keep track of seen interaction, either simple or twisted, to avoid double printing
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
                seenHashKeys.add(baseHashKey);
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


            Integer simpleNum;
            Integer twistedNum;
            if(relOriTag.equals("S")) {
                simpleNum = interaction_counts_map.get(hashKey);
                // check if there is also a twisted interaction
                String twistedHaskKey = baseHashKey.concat(";T");
                if(interaction_counts_map.containsKey(twistedHaskKey)) {
                    twistedNum = interaction_counts_map.get(twistedHaskKey);
                } else {
                    twistedNum = 0;
                }
                printStream.print("\t");
                printStream.print(simpleNum);
                printStream.print(":");
                printStream.print(twistedNum);

                if(simpleNum+twistedNum==1) {n_singleton_interactions++;}

            } else {
                twistedNum = interaction_counts_map.get(hashKey);
                // check if there is also a twisted interaction
                String simpleHaskKey = baseHashKey.concat(";S");
                if(interaction_counts_map.containsKey(simpleHaskKey)) {
                    simpleNum = interaction_counts_map.get(simpleHaskKey);
                } else {
                    simpleNum = 0;
                }
                printStream.print("\t");
                printStream.print(simpleNum);
                printStream.print(":");
                printStream.print(twistedNum);

                if(simpleNum+twistedNum==1) {n_singleton_interactions++;}

                //printStream.print(interaction_counts_map.get(hashKey).get(i));
            }
            printStream.print("\n");
            if(seenHashKeys.size()%1000000==0) {logger.trace("seenHashKeys.size(): " + seenHashKeys.size());}
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
            if(cnt%10000==0) {logger.trace(cnt);}
            Map.Entry pair = (Map.Entry)it.next();
            String hashKey = pair.getKey().toString();
            String[] fragKey = hashKey.split(";");

            // check if fragment keys exists and create array if necessary
            if(!fragment_interaction_counts_map.containsKey(fragKey[0])) {
                fragment_interaction_counts_map.put(fragKey[0], 0);
            }
            if(!fragment_interaction_counts_map.containsKey(fragKey[1])) {
                fragment_interaction_counts_map.put(fragKey[1], 0);
            }

            // increment counts for both fragments
                if(fragment_interaction_counts_map.get(fragKey[0])==0) {
                    interacting_fragment_count++;
                    if(fragmentKeyIsActive(fragKey[0])) {
                        active_interacting_fragment_count++;
                    }
                }
                Integer newVal = fragment_interaction_counts_map.get(fragKey[0]) + interaction_counts_map.get(hashKey);
                fragment_interaction_counts_map.put(fragKey[0], newVal);
                if(fragment_interaction_counts_map.get(fragKey[1])==0) {
                    interacting_fragment_count++;
                    if(fragmentKeyIsActive(fragKey[1])) {
                        active_interacting_fragment_count++;
                    }
                }
                newVal = fragment_interaction_counts_map.get(fragKey[1]) + interaction_counts_map.get(hashKey);
                fragment_interaction_counts_map.put(fragKey[1], newVal);
        }
    }


    /**
     * This function writes one tab delimited text file format to disk.
     * Each row corresponds to a restriction fragment that contains at least
     * one read of a valid pair.
     *
     * The first three fields contain the coordinates of the interacting fragment.
     * All following fields contain the number of reads (that are part of a valid pair)
     * at the interacting fragments.
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
            printStream.print("\t" + fragment_interaction_counts_map.get(hashKey));
            printStream.print("\n");
        }
    }


    /**
     * This method is primarily for testing and therefore it's public.
     *
     * @param hashKey
     *
     * @return Current number of interactions for given fragment pair.
     */
    public int getNumberOfInteractionsForKey(String hashKey) {
        int interNum = interaction_counts_map.get(hashKey);
        return interNum;
    }

    /**
     * This method is primarily for testing and therefore it's public.
     *
     * @param hashKey
     *
     * @return Current number of reads for a given fragment.
     */
    public int getNumberOfReadsAtInteractingFragmentForKey(String hashKey) {
        int interNum = fragment_interaction_counts_map.get(hashKey);
        return interNum;
    }
}
