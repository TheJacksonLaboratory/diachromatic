package org.jax.diachromatic.count;


import htsjdk.samtools.*;
import htsjdk.samtools.util.Log;

import org.jax.diachromatic.align.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *  // https://diachromatic.readthedocs.io/en/latest/mapping.html
 *
 * @author <a href="mailto:peter.robinson@jax.org">Peter Robinson</a>
 * @author <a href="mailto:peter.hansen@charite.de">Peter Hansen</a>
 * @version 0.0.0 (2018-05-03)
 *
 */
public class BadReadsCounter {
    private static final Logger logger = LoggerFactory.getLogger(Counter.class);

    /**
     * HashMap that stores interaction counts. Key: reference of digest pair; Value: SimpleTwistedCount objects.
     */
    private final Map<DigestPair,BadReadCount> dp2countsMap;

    /**
     * Stores interaction counts.
     */
    private final DigestMap digestMap;
    /**
     * A reader for the unique valid read pairs.
     */
    final private SamReader reader;

    /**
     * Iterator over reads from {@link #reader}.
     */
    final private Iterator<SAMRecord> it;

    /**
     * Keeps track of the number of processed read pairs
     */
    private int n_pairs_total = 0;

    private int n_unpaired = 0;

    private int n_same_digest = 0;

    private final Map<Digest, BadReadCount> badReadCountMap;


    /**
     Un-ligated due to size (Tag: UL)
     Un-ligated due to same digest (Tag: ULSI)
     Self-ligated due to size (Tag: SL)
     Self-ligated due to same digest (Tag: SLSI)
     Too short chimeric (Tag: TS)
     Too long chimeric (Tag: TL)
     Valid pair (Tag: VP)
     */
    private final Map<String, Integer> readPairCategoryCounts;



    public BadReadsCounter(SamReader samReader, DigestMap digestMap) {
        this.reader = samReader;
        this.digestMap = digestMap;
        this.it = reader.iterator();
        this.dp2countsMap = new HashMap<>();
        readPairCategoryCounts = new HashMap<>();
        readPairCategoryCounts.put("UL", 0);
        readPairCategoryCounts.put("ULSI", 0);
        readPairCategoryCounts.put("SI", 0);
        readPairCategoryCounts.put("SLSI", 0);
        readPairCategoryCounts.put("TS", 0);
        readPairCategoryCounts.put("TL", 0);
        readPairCategoryCounts.put("VP", 0);
        badReadCountMap = new HashMap<>();
    }

    public void countInteractions() {

        // iterate over unique valid pairs
        n_pairs_total = 0;
        while (it.hasNext()) {
            SAMRecord record1 = it.next();
            SAMRecord record2 = it.next();
            ReadPair readPair = new ReadPair(record1, record2, digestMap);
            n_pairs_total ++;
            //readPair.setRelativeOrientationTag();
            if (! readPair.isPaired()) {
                // at least one read is not uniquely mapped, skip this readpair
                n_unpaired++;
                continue;
            }
            DigestPair dp = readPair.getDigestPair();

            // String categoryTag = readPair.getCategoryTag();
            // is seemingly always NA
            String read1yy = record1.hasAttribute("YY") ? record1.getAttribute("YY").toString() : "NA";
            String read2yy = record2.hasAttribute("YY") ? record2.getAttribute("YY").toString() : "NA";
            if (! read1yy.equals(read2yy)) {
                System.out.println(read1yy + " - " + read2yy);
            } else {
                readPairCategoryCounts.merge(read1yy, 1, Integer::sum);
            }
            if (dp.isSameDigest()) {
                n_same_digest++;
                continue;
            }
            if (read1yy.equals("TS") || read1yy.equals("UL")) {
                // get too-short (TS) or un-ligated (UL) pairs only
                Digest d5 = dp.get5digest();
                Digest d3 = dp.get3digest();
                this.badReadCountMap.putIfAbsent(d5, new BadReadCount());
                this.badReadCountMap.putIfAbsent(d3, new BadReadCount());
                if (read1yy.equals("TS")) {
                    this.badReadCountMap.get(d5).incrementTs3();
                    this.badReadCountMap.get(d3).incrementTs5();
                } else {
                    this.badReadCountMap.get(d5).incrementUl3();
                    this.badReadCountMap.get(d3).incrementUl5();
                }
            }


            if (n_pairs_total % 10000000 == 0) {
                logger.trace("Number of read pairs: " + n_pairs_total);
                System.out.printf("Number of read pairs: %d\n", n_pairs_total);
                for (var e : readPairCategoryCounts.entrySet()) {
                    System.out.printf("%s: %d\n", e.getKey(), e.getValue());
                }
            }
        }
    }


    public void writeToFile(String outputName) {
        logger.info("Writing results to {}", outputName);
        System.out.printf("Writing results to %s\n", outputName);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputName))) {
            String [] headerfields = {"Digest", "5'UL", "5'TS", "3'UL", "#'TS"};
            bw.write(String.join("\t", headerfields) + "\n");
            for (var e : this.badReadCountMap.entrySet()) {
                String digest = e.getKey().toString();
                BadReadCount brc = e.getValue();
                String l = String.format("%s\t%d\t%d\t%d\t%d\n",
                        digest, brc.ul5count, brc.ts5count, brc.ul3count, brc.ts3count);
                bw.write(l );
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

