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
The default name for the output file is hicupCloneDigest.txt, but can be overridden with 
the -o option. We will use this file further below for the map command.

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
 map -b /usr/bin/bowtie2 -i /home/robinp/bin/bowtie2/hg19 -q1 seq1.fastq -q2 seq2.fastq -d hicupCloneDigest.txt  
 ```
The mapped combines mapping and filtering. The steps are
1. Call bowtie2 and map each of the input files separately (as it they were single-end sequences).
2. Then input both of the resulting SAM files. The files are filtered and
mapped accoding to the following criteria
* Discard a read pair if one or more of the reads could not be mapped by bowtei2
* Discard the pair if one of more of the reads was multimapped (bowtie2: XS:i tag exists)
* to do -- filter the read to remove typical HiC artefacts
* todo -- demand that at least one read maps to one of the VPV target regions


