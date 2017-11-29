package org.jax.diachromatic.map;


import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileWriterFactory;


import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jax.diachromatic.util.Pair;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;


public class SAMPairer {
    private static final Logger logger = LogManager.getLogger();
    String VERSION="1.0";

    private static final Log log = Log.getInstance(SAMPairer.class);

    String samPath1,samPath2;

    /** Number of unmapped forward reads (SAM flag==4) */
    int n_unmapped_read1=0;
    /** Number of unmapped reverse reads (Note: SAM flag is also 4, because the reverse reads are mapped as single-end reads). */
    int n_unmapped_read2=0;
    /** Number of pairs with 1 or 2 unmapped reads (these pairs are discarded from further analysis).*/
    int n_unmapped_pair=0;
    /** Number of forward reads that were multimapped (had an XS tag) */
    int n_multimapped_read1=0;
    /** Number of reverse reads that were multimapped (had an XS tag) */
    int n_multimapped_read2=0;
    /** Number of  read pairs where one or two reads were multimapped (had an XS tag) */
    int n_multimappedPair=0;





    public SAMPairer(String sam1, String sam2) {
        samPath1=sam1;
        samPath2=sam2;
    }

    Pair<SAMRecord,SAMRecord> getNextPair(Iterator<SAMRecord> it1,Iterator<SAMRecord> it2) {
        SAMRecord record1=null;
        SAMRecord record2=null;
        if (it1.hasNext() && it2.hasNext()) {
            record1=it1.next();
            record2=it2.next();
            return new Pair<>(record1,record2);
        } else {
            return null;
        }

    }

    public void pair() throws IOException {
        final SamReader reader1 = SamReaderFactory.makeDefault().open(new File(samPath1));
        final SamReader reader2 = SamReaderFactory.makeDefault().open(new File(samPath2));


        SAMFileHeader header = reader1.getFileHeader();
        // TODO Add CL:"/usr/bin/bowtie2-align-s --wrapper basic-0 --very-sensitive -x /home/robinp/bin/bowtie2/hg19 -p 1 - --passthrough"
        String programGroupId="@PG\tID:Diachromatic\tPN:Diachromatic\tVN:" + VERSION;//"\@PG\tID:HiCUP Mapper\tVN:" . "$hicup_module::VERSION\n";
        SAMProgramRecord programRecord = new SAMProgramRecord(programGroupId);
        header.addProgramRecord(programRecord);

        // init BAM outfile
        String outfile="test.bam";
        boolean presorted=false;
        final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header,presorted,new File(outfile));

        final ProgressLogger pl = new ProgressLogger(log, 1000000);


        Iterator<SAMRecord> it1 =reader1.iterator();
        Iterator<SAMRecord> it2 =reader2.iterator();



        Pair<SAMRecord,SAMRecord> pair = getNextPair(it1,it2);
        while (pair != null) {
            // first check whether both reads were mapped.
            int flag1 = pair.first.getFlags();
            int flag2 = pair.second.getFlags();
            //read 1 could not be aligned
            if (pair.first.getReadUnmappedFlag()) {
                n_unmapped_read1++;
                n_unmapped_pair++;
                if (pair.second.getReadUnmappedFlag()) n_unmapped_read2++;
            } else if (pair.first.getReadUnmappedFlag()) {
                n_unmapped_read2++; // note read1 must be OK if we get here...
                n_unmapped_pair++;
            } else  if ( pair.first.getAttribute("XS") != null ) {
                // Now look for multimapped reads.
                // If a read has an XS attribute, then bowtie2 multi-mapped it.
                n_multimapped_read1++;
                if (pair.second.getAttribute("XS") != null) { n_multimapped_read2++; }
                n_multimappedPair++;
            } else if (pair.second.getAttribute("XS") != null) {
                n_multimapped_read2++; // note if we are here, read1 was not multimapped
                n_multimappedPair++;
            } else {
                // If we get here, then we want to figure out where the reads map to.
                String chrom1 = pair.first.getReferenceName();
                int start1 = pair.first.getAlignmentStart();
                int end1 = pair.first.getAlignmentEnd();
                String chrom2 = pair.second.getReferenceName();
                int start2 = pair.second.getAlignmentStart();
                int end2 = pair.second.getAlignmentEnd();

                if (is_valid(chrom1, start1, end1, chrom2, start2, end2)) {
                    // do something here to add this VALID INTERACTION PAIR to a data structure
                    // also write it out to our BAM file
                    // Note we need to add corresponding bits to the SAM flag
                    pair.first.setFirstOfPairFlag(true);
                    pair.second.setSecondOfPairFlag(true);
                    System.out.println("*****  READ 1  *****");
                    SamBitflagFilter.debugDisplayBitflag(flag1);
                    System.out.println("*****  READ 2  *****");
                    SamBitflagFilter.debugDisplayBitflag(flag1);
                    // Now set the flag to indicate it is paired end data
                    System.out.println("     READ 1  read paired flag");
                    pair.first.setReadPairedFlag(true);// 0x1
                    pair.first.setProperPairFlag(true);//0x2
                    pair.second.setReadPairedFlag(true);
                    pair.second.setProperPairFlag(true);
                    SamBitflagFilter.debugDisplayBitflag(pair.first.getFlags());

                    // Indicate if pair is on the reverse strand
                    System.out.println("   READ 1&2  set mate negative strand flag");
                    pair.first.setMateNegativeStrandFlag(pair.second.getReadNegativeStrandFlag());
                    pair.second.setMateNegativeStrandFlag(pair.first.getReadNegativeStrandFlag());


                    // Set which reads are which in the pair
                    pair.first.setFirstOfPairFlag(true);
                    System.out.println("   READ 1  set first segment in template");
                    SamBitflagFilter.debugDisplayBitflag(pair.first.getFlags());
                    pair.second.setSecondOfPairFlag(true);
                    System.out.println("   READ 2  set second segment in template");
                    SamBitflagFilter.debugDisplayBitflag(pair.second.getFlags());



                    // Set the RNEXT and PNEXT values
                    // If the reference indices are the same, then the following should print "=" TODO CHeck this.
                    pair.first.setMateReferenceIndex(pair.second.getReferenceIndex());
                    pair.second.setMateReferenceIndex(pair.first.getReferenceIndex());

                    pair.first.setMateAlignmentStart(pair.second.getAlignmentStart());
                    pair.second.setMateAlignmentStart(pair.first.getAlignmentStart());

                    System.out.println(pair.first.getSAMString());

                    writer.addAlignment(pair.first);
                    writer.addAlignment(pair.second);
                }
            }


            pair=getNextPair(it1,it2);
        }
        writer.close();
    }

    /**
     * Decide if this candidate pair is valid according to the rules for capture Hi-C
     * TODO right now we only print out the pair but we need to do the right thing here!
     * @param chrom1
     * @param start1
     * @param end1
     * @param chrom2
     * @param start2
     * @param end2
     * @return
     */
    private boolean is_valid(String chrom1,int start1, int end1,String chrom2,int start2,int end2) {
        logger.trace(String.format("read 1: %s:%d-%d; read 2: %s:%d-%d",chrom1,start1,end1,chrom2,start2,end2));
        return true; // dummy function but we will do the right thing here.
    }





}
