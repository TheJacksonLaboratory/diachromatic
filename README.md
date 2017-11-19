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
