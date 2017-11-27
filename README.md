# diachromatic
Differential analysis of chromatin interactions by capture

## Summary
This package will implement parts of the capture Hi-C pipeline as done in
HiCUP, but will do so in a way that is easily integrated with VPV. The final
part of the pipeline will perform differential analysis of 
promoter reads using a likelihood ratio test.

### Preparations
For now, we will rely on bowtie2 to perform alignment. This program
can be installed on debian systems as
```
$ sudo apt-get install bowtie2
```
Currently, this installs version 2.2.6. Alternatively, download 
and install the latest version from 
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml


### Digest
The first part of the pipeline creates a digest file representing
the in silico digestion of the genome using the same enzyme(s) that
the user chose in VPV for the capture panel.

```
$ java -jar Diachromatic.jar digest -g genomeDir [-o filename]
```
When run in the digest mode, Diachromatic needs to be passed the
path to the directory that contains the genome FASTA files. This
should be the same directory used to choose the viewpoints with VPV. We will assume
that the genome has be extracted, i.e., we will look for FASTA files
in the directory and disregard gz files.
The default name for the output file is XXX, but can be overridden with 
the -o option.

We will choose to use the same output format as HiCUP, which is as follows
```
Genome:myTest   Restriction_Enzyme1:HindIII [A^AGCTT]   Restriction_Enzyme2:None        Hicup digester version 0.5.10
Chromosome      Fragment_Start_Position Fragment_End_Position   Fragment_Number RE1_Fragment_Number     5'_Restriction_Site     3'_Restriction_Site
chr12_KI270835v1_alt    1       1235    1       1       None    Re1
chr12_KI270835v1_alt    1236    2327    2       2       Re1     Re1
(...)
```
### Truncate
The next part of the pipeline truncates part of reads that contain the ligation
sequence. ToDo -- document me.

### Mapping
We now align the truncated FASTQ files to the genome. We will start for now with 
the test files from the HiCUP distribution, which are human sequences--we will 
align to the hg37 genome with bowtie 2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

We download the prebuilt index for the human hg37 genome, H. sapiens, UCSC hg19 (3.5 GB).
This needs to be unzipped. 

```$xslt
$ unzip hg19.zip 
Archive:  hg19.zip
  inflating: hg19.1.bt2              
  inflating: hg19.2.bt2              
  inflating: hg19.3.bt2              
  inflating: hg19.4.bt2              
  inflating: hg19.rev.1.bt2          
  inflating: hg19.rev.2.bt2          
  inflating: make_hg19.sh           
```
We will call the path to the directory where the index was unpacked /path/to/bowtie2index/


bowtie2 was installed using the Ubuntu apt system.

For now, we will run bowtie from the command line
```$xslt
$ bowtie2 --very-sensitive  -x /path/to/bowtie2index/hg19 --no-unal -p 1 -1 ex1.fastq -2 ex2.fastq -S outname
```



These are the steps of hicup_mapper
```
map_file....with bowtie2
$oid = open3( \*WRITER, \*READER, \*ERROR, "$config{bowtie2} --very-sensitive  -x $config{index} --no-unal -p 1 - -S $outputfile" );
open3() spawns the given $cmd and connects CHLD_OUT for reading from the child, CHLD_IN for writing to the child, and CHLD_ERR for errors.

#Write the results from the mapping into a hash
my %summary_results = extract_mapping_results("$summaryfileTemp");    #Stores the results to write to the summary file

#Now pair the files
print "Pairing files with HiCUP Mapper v$hicup_module::VERSION\n" unless $config{quiet};

open( SUMMARY, ">$summaryfile" ) or die "Could not write to $summaryfile : $!";
print SUMMARY
"File\tTotal_reads_processed\tReads_too_short_to_map\t%Reads_too_short_to_map\tUnique_alignments\t%Unique_alignments\tMultiple_alignments\t%Multiple_alignments\tFailed_to_align\t%failed_to_align\tPaired\t%Paired\n";

#Process the jobs as separate child processes (in accordance with $config{threads})
%children = ();

foreach my $fileforward ( keys %files ) {
    my $filereverse = $files{$fileforward};

    my $pid = fork();
    die "cannot fork" unless defined $pid;

    if ( $pid == 0 ) {
        pair( $fileforward, $filereverse );
        exit(0);
    } else {
        $children{$pid} = 1;
        while ( keys(%children) == $config{threads} ) {
            sleep(1);
            reaper();
        }
    }
}

#Make sure all child processes have terminated before exiting
do {
    sleep(1);
    reaper();
} until ( keys(%children) == 0 );

print "Pairing complete\n" unless $config{quiet};

#Remove intermediate *.map* files
my @map_files = fileNamer(\@filenames, \%config, 'mapper', 0, 0, 0, 1, 0);
foreach my $mapFile (@map_files) {
    unlink $config{outdir}.$mapFile or warn "Could not delete '$config{outdir}.$mapFile'\n";
}

close SUMMARY or die "Could not close filehandle on '$summaryfile' : $!";

#Produce summary graph
unless ( $config{r} eq '0' ) {    #R not installed/found
    my $command = $config{r} . 'script ' . "$Bin/r_scripts/hicup_mapper_summary.r $config{outdir} $summaryfile";
    !system("$command") or warn "Could not produce hicup_mapper summary bar chart: $command: $!";
}

exit(0);

```