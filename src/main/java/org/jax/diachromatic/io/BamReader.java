package org.jax.diachromatic.io;


import htsjdk.samtools.DefaultSAMRecordFactory;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.SeekableStream;

import java.io.File;

/**
 * Created by hansep on 11/26/17.
 */
public class BamReader {

    final private  String bamPath;

    public BamReader(String path) {
        bamPath = null;
    }

    public void openReadAndWriteBAMFile(SamReader reader) {
        for (final SAMRecord samRecord : reader) {
            System.out.println(samRecord.getReferenceName());
        }
    }




}
