package org.jax.diachromatic.map;


import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileWriterFactory;


import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.jax.diachromatic.util.Pair;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Iterator;


public class SAMPairer {

    String VERSION="1.0";

    private static final Log log = Log.getInstance(SAMPairer.class);

    String samPath1,samPath2;

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
            //
            int flag1 = pair.first.getFlags();
            int flag2 = pair.second.getFlags();

            System.out.println("flags 1. "+flag1 + ", 2. "+flag2);
        writer.addAlignment(pair.first);
        writer.addAlignment(pair.second);

            pair=getNextPair(it1,it2);
        }
        writer.close();
    }

}
