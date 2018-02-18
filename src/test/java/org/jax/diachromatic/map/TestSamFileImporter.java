package org.jax.diachromatic.map;

import org.jax.diachromatic.util.Pair;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

class TestSamFileImporter {

    private final String forwardSamPath;
    private final String reverseSamPath;
    /** A counter that helps us to number the digests. */
    private  int digestcounter=0;
    private  String restrictionsite=null;
    private  Map<String,List<Digest>> digestmap;
    private  static SAMPairer sampairer=null;
    /** We will; use this to keep track of the read pairs from the test files forwardtest and reverse test.
     * Note that we have labeled the read pairs in those file to make them easy find for these tests.
     */
    private static Map<String,ReadPair> readpairmap=null;

    private TestSamFileImporter() {
        forwardSamPath=Paths.get("src","test","resources","data","sam","forwardtest.sam").toAbsolutePath().toString();
        reverseSamPath = Paths.get("src","test","resources","data","sam","reversetest.sam").toAbsolutePath().toString();
    }

    /**
     * @return a singleton {@link SAMPairer} object with data from forward.sam and reverse.sam
     */
    static Map<String,ReadPair> retrieveReadPairsForTesting() {
        if (sampairer==null || readpairmap==null) {
            initData();
        }
        return readpairmap;
    }


   static SAMPairer retrieveSAMPairerForTesting() {
        if (sampairer==null || readpairmap==null) {
            initData();
        }
        return sampairer;
    }


    private static void initData() {
        TestSamFileImporter importer=new TestSamFileImporter();
        importer.inputFile();
        ReadPair pair;
        while ((pair = sampairer.getNextPair())!=null) {
            readpairmap.put(pair.forward().getReadName(),pair);
            System.err.println(pair.forward().getReadName());
        }
    }





    private void inputFile() {
        // make fake digests for the reads we will test
        restrictionsite="HindIII";
        digestmap=new HashMap<>();
        List<Digest> digests = makeFakeDigestList("chr1",new Pair<>(221612607,221618800));
        digestmap.put("chr1",digests);
        // chr7 -- for SRR071233.24017.dangling_ends
        digests=makeFakeDigestList("chr7",new Pair<>(34161770,34167291)	);
        digestmap.put("chr7",digests);
        // Note that we put in three contiguous segments, here. This is important for
        // test testContiguous(); if we leave out the middle segment, then the sixth read would
        // mistakenly be classified as religated (the reads are on the two outer Digests of a triplet!).
        digests = makeFakeDigestList("chr10", new Pair<>(53_177_087,53_179_904),
                new Pair<>(53_179_905,	53_180_587), new Pair<>(53_180_588,53187288));
        digestmap.put("chr10",digests);

        digests=makeFakeDigestList("chr11",new Pair<>(92314037,92316529),new Pair<>(92316530,92317574));
        digestmap.put("chr11",digests);
        // Note -- the following command will give these two digests neighboring digest numbers
        // this is needed for the test and resembles what we would have for a full dataset.
        digests=makeFakeDigestList("chr13",new Pair<>(31_421_583,31_425_191),new Pair<>(31425192,	31425873	));

        digestmap.put("chr13",digests);
        //chr15	91492555	91497580
        digests=makeFakeDigestList("chr15",new Pair<>(91_492_555,91_497_580));
        digestmap.put("chr15",digests);//chr15	91492809
        digests=makeFakeDigestList("chr16",new Pair<>(31497401,31527040), new Pair<>(84172259,84175274));
        digestmap.put("chr16",digests);
        digests=makeFakeDigestList("chr17",new Pair<>(22262265,22262874));
        digestmap.put("chr17",digests);
        digests=makeFakeDigestList("chr18",new Pair<>(71910154,71919237));
        digestmap.put("chr18",digests);
        readpairmap = new HashMap<>();
        boolean outputRejectedReads=false;
        sampairer = new SAMPairer(forwardSamPath,reverseSamPath,digestmap,outputRejectedReads);
    }


    /**
     * We require a digest list for the {@link SAMPairer} constructor. We make a "fake" digest list for
     * simplicity that will allow us to perform testing of various functionalities.
     * @param chrom chromosome
     * @param pos_pair pairs of integers with start and end position of the Digests
     * @return list of "fake" Digest object
     */
    @SafeVarargs
    private final List<Digest> makeFakeDigestList(String chrom,Pair<Integer,Integer> ...pos_pair) {
        List<Digest> dlist = new ArrayList<>();
        for (Pair<Integer,Integer> p : pos_pair) {
            Digest d = makeFakeDigest(chrom,p.first,p.second);
            dlist.add(d);
        }
        return dlist;
    }

    /**
     *
     * @param chr
     * @param topos
     * @param frompos
     * @return one fake {@link Digest} object (See {@link #makeFakeDigestList(String, Pair[])}).
     */
    private Digest makeFakeDigest(String chr, int frompos, int topos) {
        String fields[]=new String[6];
        fields[0]=chr;
        fields[1]=String.valueOf(frompos);
        fields[2]=String.valueOf(topos);
        fields[3]=String.valueOf(++digestcounter);
        fields[4]=restrictionsite;
        fields[5]=restrictionsite;
        return new Digest(fields);
    }







}
