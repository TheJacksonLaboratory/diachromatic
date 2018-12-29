package org.jax.diachromatic.allelespec;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

public class AlleleSpecificityChecker {

    /** Path to the BAM file from a capture Hi-C experiment. */
    private final String bamFilePath;
    /** SAM Reader corresponding to {@link #bamFilePath}. */
    private final SamReader samReader;


    /**
     * For the initial part of development, we will not use the digest file.
     * @param bamPath path to a BAM file that represents a CHC experiment
     * @param digestFilePath TODO
     */
    public AlleleSpecificityChecker(String bamPath,String digestFilePath){
        this.bamFilePath=bamPath;
        SamReaderFactory srf=SamReaderFactory.make();
        srf.validationStringency(ValidationStringency.LENIENT);
        this.samReader= srf.open(new File(bamPath));
    }


    public void readAlleles() {
        SAMFileHeader header = samReader.getFileHeader();
        IntervalList il = new IntervalList(header);
        String chrom= "chrY";
        int from= 1;
        int to= 100;
        //debugPrintHeader(header);
        List<SamRecordFilter> filter = ImmutableList.of(); // ?? TODO how to do this and do we want filters at all?
        SamLocusIterator samLocIter= new SamLocusIterator(samReader, il, true);
        samLocIter = new SamLocusIterator(samReader);
        int count=0;
        int n_empty=0;
        while (samLocIter.hasNext()) {
            SamLocusIterator.LocusInfo locus = samLocIter.next();
            if (locus.isEmpty()) {
                n_empty++;
                continue; // no reads at this position
            }
            //if (locus.getDeletedInRecord().size()<1) continue;
            debugPrintLocus(locus);
            List<SamLocusIterator.RecordAndOffset> recordAndOffsetList = locus.getRecordAndOffsets();
            for (SamLocusIterator.RecordAndOffset rao : recordAndOffsetList) {
               // System.out.println(rao.getBaseQualities());
            }
            if (count++>100)break;
        }
        System.out.println("We found " + n_empty + " empty positions");
    }


    private String getBaseDistribution(List<SamLocusIterator.RecordAndOffset> recordAndOffsetList){
        int A_count=0,C_count=0,G_count=0,T_count=0;
        final int QUALITY_THRESHOLD=10;
        final byte A_BASE=65;
        final byte C_BASE=67;
        final byte G_BASE=71;
        final byte T_BASE=84;
        for (SamLocusIterator.RecordAndOffset rao : recordAndOffsetList) {
            if (rao.getBaseQuality() > QUALITY_THRESHOLD) {

                switch (rao.getReadBase()) {
                    case A_BASE: A_count++;break;
                    case C_BASE: C_count++;break;
                    case G_BASE: G_count++;break;
                    case T_BASE: T_count++;break;
                    default:
                        System.out.println(rao.getReadBase() + "="+(char)rao.getReadBase());
                        System.exit(1);
                }
            }
        }
        return String.format("A:%d C:%d G:%d T:%d",A_count,C_count,G_count,T_count);
    }

    private void debugPrintLocus(SamLocusIterator.LocusInfo locus) {
        System.out.println("locus:");
        System.out.println("\tname: "+locus.getSequenceName());
        System.out.println("\tpos: "+locus.getPosition());
        System.out.println("\tsize: "+locus.size());
        System.out.println("\tisEmpty?: "+locus.isEmpty());
        System.out.println("\tsequence length: "+locus.getSequenceLength());
        System.out.println("\tsequence index: "+locus.getSequenceIndex());
        System.out.println("\ttostring: "+locus.toString());
        List<SamLocusIterator.RecordAndOffset> recordAndOffsetList =locus.getRecordAndOffsets();
        System.out.println("\trecords and offsets: ");
        System.out.println(getBaseDistribution(recordAndOffsetList));
        List<SamLocusIterator.RecordAndOffset> deleted=locus.getDeletedInRecord();
        debugPrintDeleted(deleted);
        return;/*
        for (SamLocusIterator.RecordAndOffset rao : recordAndOffsetList) {
            System.out.println("\treadname: "+rao.getReadName());
            System.out.println("\t\ttostring: "+rao.toString());
            System.out.println("\t\tlen: "+rao.getLength());
            System.out.println("\t\tbasequal: "+rao.getBaseQuality());
            System.out.println("\t\treadbase: "+(char)rao.getReadBase());
        }
    */
    }



    private void debugPrintDeleted(List<SamLocusIterator.RecordAndOffset> recordAndOffsetList){
        int A_count=0,C_count=0,G_count=0,T_count=0;
        final int QUALITY_THRESHOLD=10;
        final byte A_BASE=65;
        final byte C_BASE=67;
        final byte G_BASE=71;
        final byte T_BASE=84;
        System.out.println("Deleted n=" + recordAndOffsetList.size());
        for (SamLocusIterator.RecordAndOffset rao : recordAndOffsetList) {
            System.out.println("\tlen="+rao.getLength());
            if (rao.getBaseQuality() > QUALITY_THRESHOLD) {

                switch (rao.getReadBase()) {
                    case A_BASE: A_count++;break;
                    case C_BASE: C_count++;break;
                    case G_BASE: G_count++;break;
                    case T_BASE: T_count++;break;
                    default:
                        System.out.println(rao.getReadBase() + "="+(char)rao.getReadBase());
                        System.exit(1);
                }
            }
        }
        //return String.format("A:%d C:%d G:%d T:%d",A_count,C_count,G_count,T_count);
    }

    private void debugPrintHeader(SAMFileHeader header) {
        // printout entire SAM File header section.
        System.out.println(header.getSAMString());
    }


}
