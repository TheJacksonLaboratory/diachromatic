
.. HiCUP
.. ~~~~~
..
.. HiCUp is currently a standard tool for capture Hi-C Q/C and preprocessing and much of the preprocessing
.. for diachromatic is based on Hicup (and we cite it). Here is how Hicup was used to generate results from
.. the Hicup test dataset.
..
..
.. bowtie and index
.. ----------------
..
.. The test data set was digested with HindIII. We will align it to Hg19. We will first create an index file for
.. bowtie2 (the files used to create the index and to map the reads need to be the same!).
..
.. I used the Hg19 files from UCSC, and combined them into a single file for convenience. ::
..
..   $ cat *.fa > hg19.fa
..
.. I used the bowtie2 indexer. ::
..
..     $ bowtie2-build hg19.fa HG19
..
.. This command uses the FASTA file and specifies an output suffix of ``HG19``. (Now go drink some coffee).
..
..
.. Digestion
.. ---------
..
.. This command is performed separately from the rest of the hicup pipeline. Assuming the path to the genome fasta file
.. used to create the bowtie index is ``/path/to/hg19.fa``, then use the following command: ::
..
..     $ ./hicup_digester --re1 A^AGCTT,HindIII --genome hg19 /path/to/hg19.fa
..
.. This creates a file named ``Digest_hg19_HindIII_None_21-27-22_04-01-2018.txt`` (the time of day was used to name the file).
..
..
.. Running the pipeline
.. --------------------
..
.. With this in hand, we can run the main pipeline. The settings file is (abbreviated). ::
..
..     Outdir: mytest
..     Threads: 1
..     #Suppress progress updates (0: off, 1: on)
..     Quiet:0
..     #Retain intermediate pipeline files (0: off, 1: on)
..     Keep:1
..     #Compress outputfiles (0: off, 1: on)
..     Zip:0
..     Bowtie2: /usr/bin/bowtie2
..     #Remember to include the basename of the genome indices
..     Index: /home/peter/data/ucsc/hg19/HG19
..     #Path to the genome digest file produced by hicup_digester
..     Digest: Digest_hg19_HindIII_None_21-27-22_04-01-2018.txt
..     Format: Sanger
..     #Maximum di-tag length (optional parameter)
..     Longest: 800
..     #Minimum di-tag length (optional parameter)
..     Shortest: 150
..     #FASTQ files to be analysed, placing paired files on adjacent lines
..     test_dataset/test_dataset1.fastq
..     test_dataset/test_dataset2.fastq
..
.. With this in a file called myhicup.conf, we can run the main hicup pipeline as follows. The results of the run will be put into the ``mytest`` directory (which must be created before running the script). ::
..
..      $ ./hicup -config myhicup.conf
..
.. My goal is to create two small SAM files for testing the class SAMPairer. To do this, I commented out the following lines
.. in the hicup_mapper script. ::
..
..     foreach my $mapFile (@map_files) {
..         unlink $config{outdir}.$mapFile or warn "Could not delete '$config{outdir}.$mapFile'\n";
..     }
..
.. Sure enough, the bowtie2 single-end alignments are now retained.
..
..     * test_dataset1.map.sam
..     * test_dataset2.map.sam
..
.. These can be used in conjunction with the other output files of hicup to identify read pairs that should be filtered
.. out because of mapping issues or artefacts, as well as read pairs that are ok. We can test most of the diachromatic
.. code using a small SAM file that is excerpted from these.



.. Finding digests for testing
.. ~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. Note that many of the readpair functions require a Digest object. The following script can help to find the
.. positions of the digests. These were used in the makeFakeDigest functions in the test classes. ::
..
..     #!/usr/bin/perl -w
..     use strict;
..     use IO::File;
..
..     my $fname = shift or die "need to pass digest file name";
..     my $chr = shift or die "need to pass chromosome\n";
..     my $pos= shift or die "need to pass position";
..     print "will analyse $fname, $chr, $pos\n";
..
..     my $fh=new IO::File("$fname") or die "$!";
..     while (my $line=<$fh>) {
..         my @a=split(m/\t/,$line);
..         my $chrom=$a[0];
..         #print "chrom=$chrom and chr=$chr\n";
..         next if ($chr ne $chrom);
..         my $from =$a[1];
..         my $to=$a[2];
..         if ($pos>($from-100) && $pos < ($to+100)) {
..             print $line;
..             printf("position $pos is %d nucleotides 3' to start and %d nucleotides 5' to end of digest [len=%d]\n",($pos-$from),($to-$pos),($to-$from));
..         }
..     }



.. Test class
.. ~~~~~~~~~~
.. The main tests of the logic of the Q/C code are in SAMPairerTest. There is currently one pair of sequences
.. (in forwardtest.sam and reversetest.sam) for each of the tests we perform. ::
..
.. 	SRR071233.1     67      chr16   31526917        8       40M     =       84175204        0       NAAGATACCTTGACCGCTCATCCCCTGNNTTCATGAAAGA        !##########################!!###########        AS:i:-13  XN:i:0  XM:i:8  XO:i:0  XG:i:0  NM:i:8  MD:Z:0C26A0C6G0T0C0T0T0 YT:Z:UU
.. 	SRR071233.1     131     chr16   84175204        42      40M     =       31526917        0       AGAACCCATTCACACTCCCGCCAGCAGCAGGTTCGTGCCA        @BABA@BBBBBBBB?BBBB@:?AAAB5<BAA92A=2:;77        AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:40 YT:Z:UU
..
.. The first read should be set to 67 [read paired (0x1); read mapped in proper pair (0x2);first in pair (0x40)]. The reverse read is 131 [read paired (0x1); read mapped in proper pair (0x2); second in pair (0x80)].
..
..
.. * Test mapping
..
.. The paired FASTQ files hg19_HindIII_test_data_sam_flags_1.fast1 and hg19_HindIII_test_data_sam_flags_2.fastq were
.. processed with the command: ::
..
..     $ java -jar Diachromatic.jar map -b /usr/bin/bowtie2 -i /path-to/bowtie2-index/hg19 -q hg19_HindIII_test_data_sam_flags_1.fastq -r fastq/hg19_HindIII_test_data_sam_flags_2.fastq -d hg38digest
..
.. The resulting SAM files are being used for unit testing (to simplify and robustify testing).

