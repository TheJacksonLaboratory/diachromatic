
Mapping and categorization of Hi-C paired-end reads
===================================================

Independent mapping of forward and reverse paired-end reads using bowtie2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two reads of any given valid Hi-C read pair stem from two different interacting genomic regions that can be
separated by a large number of nucleotides within the same chromosome (**cis interactions**) or even be located on
different chromosomes (**trans interactions**). For this reason, the distance between the two 5' ends of the reads can
no longer be interpreted as the *insert size*, and the truncated forward (R1) and reverse (R2) reads have to be mapped
independently.

Diachromatic executes ``bowtie2`` separately for R1 and R2 with the ``--very-sensitive`` option. Individual reads mapping
to multiple locations are typically discarded. Diachromatic provides two levels of stringency
for the definition of multi-mapped reads:
    1. **Very stringent definition:** There is no second best alignment for the given read. In this case the line in the SAM record produced by ``bowtie2`` contains no ``XS`` tag. Use Diachromatic's ``--bowtie-stringent-unique`` or ``-bsu`` option in order to use this level of stringency.
    2. **Less stringent definition:** There can be a second best alignment, but the score of the alignment (MAPQ) needs to e greater or equal than 30 and the difference of the mapping scores between the best and second best alignment must be greater or equal than 10. This definition was borrowed from HiCUP (version v0.6.0 and higher). Diachromatic uses this option by default.


Pairing of properly mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The independently mapped reads are written to two temporary SAM files, whereby the order of read records in the
truncated FASTQ files is retained by using bowtie2's option ``--reorder``. I a next step, Diachromatic iterates
simultaneously over the two SAM files. Pairs for which at least one read could not be mapped uniquely are discarded,
whereas all other pairs are futher subdivided into different categories comprising valid interaction and artefactual
read pairs arising from shortcomings of the Hi-C protocol.

Categorization of mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hi-C fragments arise from cross-linked chromatin passing through three successive experimental processing steps:
*restriction digest*, *re-ligation* and *shearing* (see illustration below). Different fragments differ with regard to their
formation history.

.. figure:: img/fragment_formation.png
    :align: center


For Hi-C, cross-linked chromatin is digested using one or more restriction enzymes,
which results in restriction fragments whose ends re-ligate thereby forming ligation junctions.
The shearing step further increases the diversity of fragments by introducing DNA breakpoints representing a second type
of ends in addition those introduced by digestion.
Fragment ends corresponding to restriction enzyme cutting sites are generally referred to as *dangling ends* because
they failed to re-ligate.

In total, three categories of fragments are distinguished within Diachromatic: **hybrid fragments** that arise from
re-ligation between ends of different restriction fragments and two artifact types that correspond to single
restriction fragments whose ends failed to re-ligate with other fragments, either because both ends remained **un-ligated**
or **self-ligated** with each other. Hybrid fragments correspond to valid interactions but also to cross-ligation
artifacts depending on whether the re-ligtion occurred within the same protein-DNA complex or between different complexes.
Paired-end sequencing of hybrid fragments may results in all possible relative orientations, i.e. reads of given pairs
may pointing *inwards*, *outwards* or in the *same direction*.
In contrast to that, sequencing of un-ligated fragments results in inward pointing pairs only, and sequencing of
self-ligated fragments results in outward pointing pairs only.

Due to the fact that the read pair orientations overlap for the different categories, the identification of read pairs
arising from un-ligated or self-ligated fragments additionally requires the definition of a size threshold that
corresponds to the **average size of fragments of the Hi-C library**.
Roughly speaking, the underlying idea is that such pairs can be assumed to span distances not much larger than this
threshold only.
However, the determination of the size of a given fragment is not straightforward for Hi-C for several reasons.
First, the size has to be calculated differently depending the category of the fragment.
For un-ligated fragments, the size corresponds to the distance between the 5’ end mapping positions of the two reads as
usual, whereas for self-ligation and hybrid fragments the size is calculated as the sum of the two distances between
the 5' ends of the mapped reads and the next occurrence of a cutting motif in 3' direction which is assumed to correspond
to the ligation junction (Wingett 2015).
The problem with this approach is that in fact the ligation junction cannot be unambiguously determined, because the
digestion of genome is not necessarily complete, i.e. there may be restriction fragments containing uncut restriction
sites (in the illustration marked with asterisk).
In such cases, the size is underestimated, because, for lack of further information, simply the first occurrence of a cutting
motif is interpreted as the one that corresponds to the ligation junction.
Therefore, Diachromatic does not use this approach but requires this parameter to be specified using the ``-l <size>`` option.
We recommend to use external tool such as the `peak caller Q`_ for fragment size estimation.
Even though the Hamming distance method implemented in Q is intended for ChIP-seq data, it is also suitable for Hi-C,
because at restriction sites, the reads distribute in a strand specific fashion that is similar to that observed for
ChIP-seq reads. Within Diachromatic, inward pointing read pairs for which the distance between the 5' ends is less than
the specified threshold are categorized as un-ligated pairs, whereas outward pointing read pairs with a calculated size
smaller than the threshold are categorized as self-ligated pairs.

.. _peak caller Q: http://charite.github.io/Q/

Valid read pairs arising from genuine chromatin-chromatin interactions cannot be distinguished from those arising from
**cross-ligation** events.
However, the overall extend of **cross-ligation** is estimated for given experiments.
Based on the assumption that cross-ligation between DNA fragments of different chromosomes (trans) occurs more likely
as compared to cross-ligation between DNA fragments of the same chromosome (cis), the ratio of the numbers of cis
and trans read pairs is taken as an indicator of poor Hi-C libraries that contain lots of false positive interaction
pairs arising from spurious cross-ligation events (Wingett 2015, Nagano 2015).
However, it has been pointed out that this quality measure depends also on other factors such as the genome size and
number of chromosomes of the analyzed species (Wingett 2015). Diachromatic provides an alternative and more robust quality metric that
can be used to access the extent of cross-ligation. Amongst the trans read pairs, we generally observe a large proportion
of restriction fragments that are connected by single read pairs only. The number of all possible different cross-ligation
events (including cis and trans) can roughly be estimated as the square number of all restriction fragments across the
entire genome. Given this huge number, we reasoned that it is very unlikely that the same cross-ligation event occurs
twice. Therefore, we defined a **cross-ligation coefficient (CLC)** as the ratio of singleton read pairs and all read pairs.

Running Diachromatic's align subcommand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the alignment and counting step. ::

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


Output files
~~~~~~~~~~~~

Two output files will be produced:

    * ``prefix.valid.bam`` contains all uniquely mapped pairs. Known artifacts and duplicated reads are removed. This file can be used for downstream analyses.

    * ``prefix.rejected.bam`` contains all pairs that show characteristics of known artifacts:

        * insert too long (Tag: ``TB``)
        * insert too short (Tag: ``TS``)
        * circularized read (Tag: ``SL``)
        * same dangling end (Tag: ``DE``)
        * same internal (Tag: ``SI``)
        * re-ligation (Tag: ``RL``)
        * contiguous (Tag: ``CT``)

    * ``prefix.align.stats.``

Read pairs for which one read cannot be mapped or cannot be mapped uniquely (bowtie2: XS:i tag exists) will be discarded completely. Statistics about the numbers of unmappable reads, multimappable reads, and artifact pairs will be written to the screen.




