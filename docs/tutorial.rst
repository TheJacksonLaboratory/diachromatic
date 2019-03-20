Tutorial
========

Here, we present a complete tutorial for using Diachromatic for processing and quality control of Capture Hi-C reads.


The test dataset
~~~~~~~~~~~~~~~~
Note: for now I am using the HiCUP test dataset, but we should create our own. I am describing the entire tutorial using the
HiCUP dataset for now TODO make this up to date!

To get the data, visit the `download site <https://www.bioinformatics.babraham.ac.uk/projects/download.html#hicup>`_ and
download the HiCUP_test_dataset. Extract it using this command. ::

    $  tar xvfz test_dataset.tar.gz

This will create a directory called ``test_dataset``, which we will symbolize with ``TESTDIR`` in the following examples.



Truncation
~~~~~~~~~~
The first step of processing with data with Diachromatic is processing of the raw FASTQ files to recognize and truncate
reads with filled-in ligation juctions, which indicate reads that include the junction of the chimeric CHC fragment. If
a single read is chimeric, it is not possible to map it to one locus, and therefore the 3' portion of the chimeric read
is removed ("truncated"), leaving behind the 5' portion of the read that should map to a specific locus. If the 5' sequence
is too short to be mapped, the entire read pair is discarded. This is performed with the truncate command.::


    $ java -jar Diachromatic.jar truncate \
        -q ${TESTDIR}/test_dataset1.fastq \
        -r ${TESTDIR}/test_dataset2.fastq \
        -e HinDIII
        -o HinD3
        -x foo

In practice, only about XXXX percent of the readpairs are truncated.


Mapping
~~~~~~~
The second step of the pipeline is to map the truncated read pairs to the target genome. Diachromatic uses bowtie2 to perform the
mapping, and then creates a BAM file containing the valid read pairs. If desired, Diachromatic also outputs BAM files
with the discarded (arterfactual or unmappable reads).

To use the mapping file, we need to have a file that shows the locations of restriction digests across the genome.
Diachromatic users can use GOPHER to create probes and the digest file TODO -- MORE DETAIL.

It is recommended that users download bowtie index files from the bowtie site. In the following, we have
downloaded and extracted these files to a directory called ``/data/bt_indices`` (in this directory, there are multiple index files
hg19.1.bt2, hg19.3.bt2, hg19.rev.1.bt2, hg19.2.bt2, hg19.4.bt2, hg19.rev.2.bt2).


Use the following command to run the alignment step. ::

    $ java -jar target/Diachromatic.jar align \
        -b /usr/bin/bowtie2 \
        -i /data/bt_indices/hg19 \
        -q foo.truncated_R1.fastq.gz \
        -r foo.truncated_R2.fq.gz \
        -d hg19_HinDIII_DigestedGenome.txt




This document contains a summary of how some of the unit tests were designed. Probably this page can/should be
deleted before publication.



Summarize
~~~~~~~~~
To run the Summarize command with the truncate data, run the following command. ::

    $ java -jar target/Diachromatic.jar summarize \
        -o HinD3 \
        -x foo \
        -t HinD3/foo.truncation.stats


This will generate an HTML file called ``HinD3/foo.summary.stats.html``.


HiCUP
~~~~~
HiCUp is currently a standard tool for capture Hi-C Q/C and preprocessing and much of the preprocuessing
for diachromatic is based on Hicup (and we cite it). Here is how Hicup was used to generate results from
the Hicup test dataset.

bowtie and index
~~~~~~~~~~~~~~~~
The test data set was digested with HindIII. We will align it to Hg19. We will first create an index file for
bowtie2 (the files used to create the index and to map the reads need to be the same!).

I used the Hg19 files from UCSC, and combined them into a single file for convenience. ::

  $ cat *.fa > hg19.fa

I used the bowtie2 indexer. ::

    $ bowtie2-build hg19.fa HG19

This command uses the FASTA file and specifies an output suffix of ``HG19``. (Now go drink some coffee).

Digestion
~~~~~~~~~
This command is performed separately from the rest of the hicup pipeline. Assuming the path to the genome fasta file
used to create the bowtie index is ``/path/to/hg19.fa', then use the following command. ::

    $ ./hicup_digester --re1 A^AGCTT,HindIII --genome hg19 /path/to/hg19.fa

This creates a file named ``Digest_hg19_HindIII_None_21-27-22_04-01-2018.txt`` (the time of day was used to name the file).


Running the pipeline
~~~~~~~~~~~~~~~~~~~~
With this in hand, we can run the main pipeline.

The settings file is (abbreviated). ::

    Outdir: mytest
    Threads: 1
    #Suppress progress updates (0: off, 1: on)
    Quiet:0
    #Retain intermediate pipeline files (0: off, 1: on)
    Keep:1
    #Compress outputfiles (0: off, 1: on)
    Zip:0
    Bowtie2: /usr/bin/bowtie2
    #Remember to include the basename of the genome indices
    Index: /home/peter/data/ucsc/hg19/HG19
    #Path to the genome digest file produced by hicup_digester
    Digest: Digest_hg19_HindIII_None_21-27-22_04-01-2018.txt
    Format: Sanger
    #Maximum di-tag length (optional parameter)
    Longest: 800
    #Minimum di-tag length (optional parameter)
    Shortest: 150
    #FASTQ files to be analysed, placing paired files on adjacent lines
    test_dataset/test_dataset1.fastq
    test_dataset/test_dataset2.fastq

With this in a file called myhicup.conf, we can run the main hicup pipeline as follows. ::

     $ ./hicup -config myhicup.conf

The results of the run will be put into the ``mytest`` directory (which must be created before running the script).


My goal is to create two small SAM files for testing the class SAMPairer. To do this, I commented out the following lines
in the hicup_mapper script. ::

    foreach my $mapFile (@map_files) {
        unlink $config{outdir}.$mapFile or warn "Could not delete '$config{outdir}.$mapFile'\n";
    }

Sure enough, the bowtie2 single-end alignments are now retained.

    * test_dataset1.map.sam
    * test_dataset2.map.sam

These can be used in conjunction with the other output files of hicup to identify read pairs that should be filtered
out because of mapping issues or artefacts, as well as read pairs that are ok. We can test most of the diachromatic
code using a small SAM file that is excerpted from these.

Finding digests for testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Note that many of the readpair functions require a Digest object. The following script can help to find the
positions of the digests, these were used in the makeFakeDigest functions in the test classes


    #!/usr/bin/perl -w
    use strict;
    use IO::File;

    my $fname = shift or die "need to pass digest file name";
    my $chr = shift or die "need to pass chromosome\n";
    my $pos= shift or die "need to pass position";
    print "will analyse $fname, $chr, $pos\n";

    my $fh=new IO::File("$fname") or die "$!";
    while (my $line=<$fh>) {
        my @a=split(m/\t/,$line);
        my $chrom=$a[0];
        #print "chrom=$chrom and chr=$chr\n";
        next if ($chr ne $chrom);
        my $from =$a[1];
        my $to=$a[2];
        if ($pos>($from-100) && $pos < ($to+100)) {
            print $line;
            printf("position $pos is %d nucleotides 3' to start and %d nucleotides 5' to end of digest [len=%d]\n",($pos-$from),($to-$pos),($to-$from));
        }
    }





Test class
~~~~~~~~~~
The main tests of the logic of the Q/C code are in SAMPairerTest. There is currently one pair of sequences
(in forwardtest.sam and reversetest.sam) for each of the tests we perform.

SRR071233.1
SRR071233.1     67      chr16   31526917        8       40M     =       84175204        0       NAAGATACCTTGACCGCTCATCCCCTGNNTTCATGAAAGA        !##########################!!###########        AS:i:-13
        XN:i:0  XM:i:8  XO:i:0  XG:i:0  NM:i:8  MD:Z:0C26A0C6G0T0C0T0T0 YT:Z:UU
SRR071233.1     131     chr16   84175204        42      40M     =       31526917        0       AGAACCCATTCACACTCCCGCCAGCAGCAGGTTCGTGCCA        @BABA@BBBBBBBB?BBBB@:?AAAB5<BAA92A=2:;77        AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU

The first read should be set to 67 [read paired (0x1); read mapped in proper pair (0x2);first in pair (0x40)]. The reverse read is
131 [read paired (0x1); read mapped in proper pair (0x2); second in pair (0x80)].

* Test mapping


The paired FASTQ files hg19_HindIII_test_data_sam_flags_1.fast1 and hg19_HindIII_test_data_sam_flags_2.fastq were
processed with the command

    $ java -jar Diachromatic.jar map -b /usr/bin/bowtie2 -i /path-to/bowtie2-index/hg19 -q hg19_HindIII_test_data_sam_flags_1.fastq -r fastq/hg19_HindIII_test_data_sam_flags_2.fastq -d hg38digest

The resulting SAM files are being used for unit testing (to simplify and robustify testing).