Mapping and categorization of Hi-C paired-end reads
===================================================

Independent mapping of forward and reverse paired-end reads using bowtie2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two reads of a valid Hi-C read pair come from two different interacting genomic regions that can be
separated by a large number of nucleotides within the same chromosome (**cis interactions**) or even be located on
different chromosomes (**trans interactions**). For this reason, the distance between the two 5' ends of the reads can
no longer be interpreted as the *insert size*, and the truncated forward (R1) and reverse (R2) reads have to be mapped
independently.

Diachromatic executes ``bowtie2`` separately for R1 and R2 with the ``--very-sensitive`` option. Individual reads mapping
to multiple locations are typically discarded. Diachromatic provides two levels of stringency
for the definition of multi-mapped reads:
    1. **Very stringent mapping:** There is no second best alignment for the given read. In this case the line in the SAM record produced by ``bowtie2`` contains no ``XS`` tag. Use Diachromatic's ``--bowtie-stringent-unique`` or ``-bsu`` option in order to use this level of stringency.
    2. **Less stringent mapping:** There can be a second best alignment, but the score of the alignment (MAPQ) needs to e greater or equal than 30 and the difference of the mapping scores between the best and second best alignment must be greater or equal than 10 (c.f. `HiCUP <https://www.bioinformatics.babraham.ac.uk/projects/hicup/>`_). Diachromatic uses this option by default.


Pairing of properly mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The independently mapped reads are written to two temporary SAM files, whereby the order of read records in the
truncated FASTQ files is retained by using bowtie2's option ``--reorder``. In the next step, Diachromatic iterates
simultaneously over the two SAM files. Pairs for which at least one read could not be mapped uniquely are discarded,
and all other pairs are futher subdivided into valid and artefactual read pairs.

Categorization of mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hi-C fragments arise from cross-linked chromatin that is processed in three successive experimental steps:
*restriction digest*, *re-ligation* and *shearing* (see illustration below). It is important to understand these steps
in order to understand how Diachromatic determines is a readpair is valid or artefactual.


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
artifacts depending on whether the re-ligation occurred within the same protein-DNA complex or between different complexes.
Paired-end sequencing of hybrid fragments may results in all possible relative orientations, i.e. reads of given pairs
may pointing *inwards*, *outwards* or in the *same direction*.
In contrast to that, sequencing of un-ligated fragments results in inward pointing pairs only, and sequencing of
self-ligated fragments results in outward pointing pairs only. Due to the fact that the read pair orientations overlap
for the different categories, the identification of read pairs arising from un-ligated or self-ligated fragments
additionally requires the definition of size thresholds.


Size and size threshold for un-ligated fragments
------------------------------------------------

Reads arising from un-ligated fragments must point inwards and the size of the sequenced fragment corresponds to
the distance between the 5' end positions of the two reads. This distance is here referred to as *d*.



We recommend to use an un-ligation threshold T\ :sub:`d` that corresponds to the **average size of fragments of the Hi-C library**.



The size distribution of un-ligated fragments should be the same as for hybrid fragments.

However, the determination of this distribution is not straightforward.

For hybrid fragments the size *d'* is calculated as the sum of the two distances between the 5' ends of the mapped reads and
the next occurrence of a cutting motif in 3' direction which is assumed to correspond to the ligation junction
(`Wingett 2015 <https://www.ncbi.nlm.nih.gov/pubmed/26835000/>`_).

The problem with this approach is that in fact the ligation junction cannot be unambiguously determined,
because the digestion of genome is not necessarily complete, i.e. there may be restriction fragments containing uncut
restriction sites (in the illustration marked with asterisk).

In such cases, the size of hybrid fragments is underestimated, because, for lack of further information, simply the
first occurrence of a cutting motif is interpreted as the one that corresponds to the ligation junction.

Therefore, Diachromatic does not use this approach but requires this parameter to be specified
using the ``-l <size>`` option.

We recommend to use an external tool such as the `peak caller Q`_ for fragment size
estimation. Even though the Hamming distance method implemented in Q is intended for ChIP-seq data, it is also suitable
for Hi-C, because at restriction sites, the reads distribute in a strand specific fashion that is similar to that
observed for ChIP-seq reads.

With Diachromatic, inward pointing read pairs for which the distance between the 5' ends *d*
is less than the specified threshold *Tu* are categorized as un-ligated pairs.


Size threshold for self-ligated fragments
-----------------------------------------

Unlike reads arising from un-ligated fragments, reads arising from self-ligated must point outwards.

Furthermore, self-ligating fragments have a different size distribution than hybrid and un-ligated fragments.

The relevant size is no longer the size of the sequenced fragments that results from sonication but the
favourable size at which fragments tend to self-ligate.

Very short fragments might not self-ligate because of steric hindrance, whereas the ends of very long fragments might
be unlikely to become located in sufficient physical proximity in order to ligate.

Within Diachromatic, the size of self-ligating fragments is calculated as the sum *ds=d+d'*,
where *d* is the distance between the 5' end positions of the two reads, and *d'* is the sum of the two distances between
the 5' ends of the mapped reads and the next occurrence of a cutting motif in 3' direction (see illustration
below).

Within Diachromatic, outward pointing read pairs for which the calculated size of the self-ligating fragment *ds* is
less than the specified threshold *Ts* are categorized as self-ligated pairs.


.. figure:: img/fragment_categories.png
    :align: center

.. _peak caller Q: http://charite.github.io/Q/


Quality metrics
~~~~~~~~~~~~~~~

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

    $ java -jar target/Diachromatic.jar map -b /usr/bin/bowtie2 -i /data/bt_indices/hg38 -q prefix.truncated_R1.fq.gz -r prefix.truncated_R2.fq.gz -d hg38_DpnII_DigestedGenome.txt


+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| Short option | Long option          | Example                                                | Required | Description                                                          | Default |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -q           | --fastq-r1           | prefix.truncated_R1.fq.gz                              | yes      | Path to the truncated forward FASTQ file.                            |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -r           | --fastq-r2           | prefix.truncated_R2.fq.gz                              | yes      | Path to the truncated reverse FASTQ file.                            |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -b           | --bowtie2            | /tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2            | yes      | Path to bowtie2 executable.                                          |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -i           | bowtie2-index        | /data/indices/bowtie2/hg38/hg38                        | yes      | Path to bowtie2 index of the corresponding genome.                   |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -d           | --digest-file        | /data/GOPHER/hg38_DpnII_DigestedGenome.txt             | yes      | Path to the digest file produced with GOPHER.                        |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -od          | --out-directory      | cd4v2                                                  | no       | Directory containing the output of the align subcommand.             | results |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -op          | ---out-prefix        | stim_rep1                                              | no       | Prefix for all generated files in output directory.                  | prefix  |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -p           | --thread-num         | 15                                                     | no       | Number of threads used by bowtie2.                                   | 1       |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -j           | --output-rejected    | --                                                     | no       | If set, a BAM file containing the reject read pairs will be created. | false   |
+--------------+----------------------+--------------------------------------------------------+----------+----------------------------------------------------------------------+---------+




Output files
~~~~~~~~~~~~

The default name of the BAM file containing all unique valid pairs that can be used for downstream analysis is:

    * ``prefix.valid_pairs.aligned.bam``


If ``--output-rejected`` is set, there will be second BAM file cointaing all rejected pairs:

    * ``prefix.rejected_pairs.aligned.bam``

The optional fields of the SAM records contain information about the reasons for rejection:

    * insert too long (Tag: ``TB``)
    * insert too short (Tag: ``TS``)
    * circularized read (Tag: ``SL``)
    * same dangling end (Tag: ``DE``)
    * same internal (Tag: ``SI``)
    * re-ligation (Tag: ``RL``)
    * contiguous (Tag: ``CT``)

Furthermore, there is an ``RO`` attribute that gives the relative orientation of the pair (``R1F2``, ``R2F1``, etc.).

In addition, a file

    * ``prefix.align.stats.``

is produced that contains summary statistics about the alignment step.
