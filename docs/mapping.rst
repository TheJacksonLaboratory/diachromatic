
Mapping of paired-end Hi-C reads
================================

Preparation of the bowtie2 index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The prebuilt ``bowtie2`` indices for human hg19 (3.5 GB) and other genome builds can be downloaded from the
`bowtei2 website`_. Move the downloaded archive to an appropriate on your computer and unpack with: ::

    $ unzip hg19.zip
        Archive:  hg19.zip
        inflating: hg19.1.bt2
        inflating: hg19.2.bt2
        inflating: hg19.3.bt2
        inflating: hg19.4.bt2
        inflating: hg19.rev.1.bt2
        inflating: hg19.rev.2.bt2
        inflating: make_hg19.sh

We will call the path to the directory where the index was unpacked **/path/to/bowtie2index/**.

.. _bowtei2 website: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml


Independent mapping of forward and reverse paired-end reads using bowtie2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two reads of any given valid Hi-C read pair stem from two different interacting genomic regions that can be
separated by a large number of nucleotides within the same chromosome (**cis interactions**) or even be located on
different chromosomes (**trans interactions**). For this reason, the distance between the two 5' ends of the reads can
no longer be interpreted as the *insert size*, and the forward (R1) and reverse (R2) reads have to be mapped
independently.

Diachromatic executes ``bowtie2`` two times with the ``--very-sensitive`` option. Individual reads mapping to multiple locations
are typically discarded. Diachromatic provides two levels of stringency
for the definition of multi-mapped reads:
    1. *Very stringent definition:* There is no second best alignment for the given read. In this case the line in the SAM file produced by ``bowtie2`` contains no ``XS`` tag. Use Diachromatic's ``--bowtie-stringent-unique`` or ``-bsu`` option in order to use this level of stringency.
    2. *Less stringent definition:* There can be a second best alignment, but the score of the alignment needs to e greater than 30 and the difference of the mapping scores between the best and second best alignment must be greater than 10. This definition was adopted from HiCUP (since v0.6.0). Diachromatic uses this option by default.


Pairing of proper mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The independently mapped reads are written to two temporary SAM files, whereby the order of reads is the same for both
files, i.e. two reads of any given line consitute a pair. I a next step Diachromatic iterates simultaneously over the
two SAM files. Only pairs for which both reads could be uniquely mapped are retained and all other pairs are discarded.

Categorization of proper mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Depending on the different formation processes, Diachromatic takes into consideration five different categories
of reads pairs:

    1. Valid interactions between regions on the same chromosome (cis)
    2. Self-ligation
    3. Cross-ligation of cross-linked protein-DNA complexes from the same chromosome
    4. Valid interactions between chromosomes (trans)
    5. Cross-ligation of cross-linked protein-DNA complexes from different chromosomes

We found no criterion that could be used in order to distinguish read pairs that emerged from cross-ligation events
from others. However, we generally notice a large fraction of trans *interactions* between pairs of restriction
fragments consisting of only one read pair. We believe that those read pairs mainly result from cross-ligation events
and use their total number in order to calculate a global cross-ligation coefficient (CLC).

We also found no accurate way to distinguish read pairs that emerged from very short valid range contacts and
self-ligation events.

Instead, we estimate the average size of chimeric fragments and use this size as a self-ligation threshold.


Performing the mapping and Q/C step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use the following command: ::

    $ java -jar Diachromatic.jar align -b <bowtie2> -i <bowtie2-index> -q <fastq1> -r <fastq2> -d <digest> [-o <outfile>]

The meaning of the options is:
    * -b <bowtie2> Path to the bowtie2 executable
    * -i <bowtie2-index> Path to the bowtie2 index for the genome used to map the FASTQ files
    * --q <fastq1> Name and path to the *truncated* "forward" FASTQ file (produced in previous step)
    * --r <fastq2> Name and path to the *truncated* "reverse" FASTQ file (produced in previous step)
    * -d <digest> Path to the digest file produced in the first step
    * [-o <outfile>] This flag is optional and if it is not passed, the default name of ``diachromatic-processed.bam`` will be used.
    * [-x] If this is option is used a set, an additional BAM file for rejected pairs will be created. The general tag for rejected reads is ``YY``. See below for tags of individual artifacts.

For instance, the following command will use bowtie2 to map the two FASTQ files of a paired-end run independently (as it they were single-end sequences). Subsequently, the two resulting mappings will be paired, and pairs that show characteristics of known artifacts will be counted and sorted out. Finally, duplicates will be removed. ::

    $ java -jar target/diachromatic-0.0.2.jar map -b /usr/bin/bowtie2 -i btindex/hg19 -q hindIIIhg19chc/test_dataset1.hindIIIhg19.fastq -r hindIIIhg19chc/test_dataset2.hindIIIhg19.fastq -d hg19HindIIIdigest.txtr -o hindIII

Two output files will be produced:
    * ``diachromatic.valid.bam`` contains all uniquely mapped pairs. Known artifacts and duplicated reads are removed. This file can be used for downstream analyses.
    * ``diachromatic.rejected.bam`` contains all pairs that show characteristics of known artifacts:
        * insert too long (Tag: ``TB``)
        * insert too short (Tag: ``TS``)
        * circularized read (Tag: ``SL``)
        * same dangling end (Tag: ``DE``)
        * same internal (Tag: ``SI``)
        * re-ligation (Tag: ``RL``)
        * contiguous (Tag: ``CT``)

Read pairs for which one read cannot be mapped or cannot be mapped uniquely (bowtie2: XS:i tag exists) will be discarded completely. Statistics about the numbers of unmappable reads, multimappable reads, and artifact pairs will be written to the screen.


todo -- demand that at least one read maps to one of the VPV target regions


