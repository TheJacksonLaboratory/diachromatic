.. _rstalign:

========================================
Mapping and categorization of Hi-C reads
========================================

The two reads of a valid Hi-C read pair come from two different interacting genomic regions that can be separated by a large number of nucleotides on the same chromosome (cis) or even be located on different chromosomes (trans). The truncated forward (R1) and reverse (R2) reads have to be mapped independently.


Independent mapping of forward and reverse paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diachromatic separately executes ``bowtie2``  with the ``--very-sensitive`` option for the truncated R1 and R2 reads. Read pairs for which at least one read cannot be mapped uniquely are discarded. Diachromatic provides two levels of stringency for the definition of multi-mapped reads:

1. **Very stringent mapping:** There is no second best alignment for the given read. In this case the line in the SAM record produced by ``bowtie2`` contains no ``XS`` tag. Use Diachromatic's ``--bowtie-stringent-unique`` or ``-bsu`` option in order to use this level of stringency.
2. **Less stringent mapping:** There can be a second best alignment, but the score of the alignment (MAPQ) must be at least 30 and the difference of the mapping scores between the best and second best alignment must be at least 10 (following the recommendation of `HiCUP <https://www.bioinformatics.babraham.ac.uk/projects/hicup/>`_). Diachromatic uses this option by default.


Pairing of properly mapped read pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The independently mapped reads are written to two temporary SAM files, whereby the order of read records in the truncated FASTQ files is retained by using bowtie2's option ``--reorder``. In the next step, Diachromatic iterates simultaneously over the two SAM files. Read pairs for which both reads can be mapped uniquely are paired, i.e. the two SAM records for single-end reads are combined into one paired-end record with appropriate SAM flags reflecting the relative orientation of the reads.


Categorization of read pairs
----------------------------

Diachromatic distinguishes several read pair categories: (A) Trans reads by definition are chimeric fragments and may represent valid biological interactions or random cross-ligation events. (B) Pairs mapping to different strands of the same chromosome may originate from un-ligated or self-ligated digests. (C) Inward pointing pairs that map to the same digest must have originated from un-ligated fragments. Size thresholds are applied to the remaining fragments to categorize them as valid or artefactual. (D) Outward pointing read pairs that map the same digest must have originated from self-ligated digests. Size thresholds are applied to the remaining fragments to categorize them as valid or artefactual. (E) Read pairs mapping to the same strand can only be chimeric. However, we observe very small proportions of read pairs that are mapped to the same strand and digest. Such read pairs are classified as strange internal.

.. figure:: img/categorization.png
    :align: center


**1. Un-ligated:** The read pair points inwards and the distance between the two 5' end positions d\ :sub:`u` is smaller than T1\ :sub:`max` or both reads map to the same digest.

**2. Self-ligated:** The read pair points outwards and the calculated size of self-ligating fragments d\ :sub:`s` is smaller than a predefined self-ligation threshold T2\ :sub:`max` (Default: 3000) or both reads map to the same digest.

**3. Short chimeric:** The read pair is not in the un-ligated or self-ligated category and the calculated size d\ :sub:`c` is smaller than a specified lower threshold threshold T1\ :sub:`min` (Default: 50).

**4. Long chimeric:** The read pair  is not in the un-ligated or self-ligated category and the calculated size d\ :sub:`c` is greater than a specified lower threshold T1\ :sub:`max` (Default: 800).

**5. Valid (chimeric):** All remaining chimeric read pairs.


The decision as to whether a read-pair is valid or not is made according to:

**1.** Read pairs that map to different chromosomes or to the same strand cannot originate from un-ligated or self-ligated fragments. Therefore, they are categorized as chimeric read pairs that are valid, if the size d\ :sub:`s` is within the specified range.

**2.** Read pairs that point inwards might originate from un-ligated fragments. In such cases, the distance between the 5' end positions of the mapped reads d\ :sub:`u` corresponds to the size of the  sequenced fragment. In order to assign read pairs to the un-ligated category, we use an upper size threshold T\ :sub:`1` that should reflect the maximum plausible size of sheared fragments. Furthermore, inward pointing read pairs that map to the same digest are categorized as un-ligated.

**3.** Read pairs that point outwards might originate from self-ligated fragments. In such cases, the size d\ :sub:`s` of the potentially underlying self-ligated fragment is calculated as described above, and compared to an upper size threshold T\ :sub:`2` for self-ligated fragments. Outward pointing read pairs with d\ :sub:`s` smaller than T\ :sub:`2` are assigned to the self-ligated category. Furthermore, outward pointing read pairs that map to the same digest are categorized as self-ligated.

**4.** Read pairs arising from chimeric fragments (not un- or self-ligated) are further distinguished. Read pairs with size d\ :sub:`s` outside the specified size range of sheared fragments will be categorizesd as too small or too large, and all remaining read pairs are categorized as valid.


Dangling end read pairs
-----------------------

Fragment ends that corresponding to restriction enzyme cutting sites are referred to as dangling ends.
In theory, fragments of all categories may have dangling ends. Therefore, there is no separate class for dangling ends.
However, the number of dangling end read pairs within each of the five disjoint categories is determined and reported.


.. Dichromatic vs. HiCUP categories
.. --------------------------------
..
.. When HiCUP is executed with the ``--keep`` flag, it will create a directory containing BAM files for the individual read pair
.. categories. We applied HiCUP to the associated test data, converted the BAM files back to FASTQ format
.. and applied Diachromatic to the FASTQ files.
..
.. The following table shows the numbers of read pairs within the categories of HiCUP and Diachromatic.
..
.. For instance, HiCUP categorized 13,760 read pairs as *same internal* and 13,722 of these are uniquely mapped using Diachromatic.
.. The small differences between these numbers may be due to different bowtie versions or settings.
.. 13,645 of these uniquely mapped read pairs are categorized as un-ligated, which is the correct category for those read pairs
.. because according to our logic *same internal* read pairs correspond to un-ligated fragments.
.. However, in total 77 *same internal* read pairs are categorized as *chimeric* read pairs, which is contradictory.
.. Further investigation revealed that the 5' end positions of those read pairs are indeed mapped to the same digest but also to the same strand.
.. According to the logic implemented in Diachromatic (see decision tree) read pairs mapped to the same strand are automatically categorized as chimeric,
.. because the concept of the Hi-C fragment formation cannot explain such read pairs.
..
.. The next HiCUP category is *re-ligation*. For Diachromatic, all 1060 read pairs are mapped uniquely.
.. 58 read pairs are categorized as *self-ligated*. Further investigation of these read pairs revealed that all pairs
.. are outward pointing, which is correct for *self-ligated* pairs.
.. The 5 *re-ligation* read pairs that are categorized as *chimeric too short* are outward pointing as well but d\ :sub:`u` is greater than
.. the self-ligation threshold. However, the calcluated size calculated d\ :sub:`c` is smaller than lower threshold for sheared fragments.
..
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **HiCUP** (rows) vs. **Diachromatic** (columns)      | **# Processed pairs** | **# Uniquely mapped pairs** | **# Un-ligated** | **# Self-ligated** | **# Chimeric too short** | **# Chimeric too long** | **# Valid** |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same internal**                                    |                13,760 |                      13,722 |       **13,645** |                  0 |                       13 |                      39 |          25 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Re-ligation**                                      |                 1,060 |                       1,060 |          **842** |                 58 |                        5 |                      49 |         106 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Contiguous**                                       |                    58 |                          58 |           **53** |                  0 |                        1 |                       0 |           4 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same circularised**                                |                   428 |                         428 |                3 |            **425** |                        0 |                       0 |           0 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Wrong size**                                       |                10,321 |                      10,267 |                2 |                  0 |                **1,003** |               **9,181** |          81 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Valid**                                            |                25,915 |                      25,851 |                1 |                  5 |                      290 |                       6 |  **25,549** |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. |                                                      |                       |                             |                  |                    |                          |                         |             |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same dangling ends**                               |                 2,475 |                       2,473 |        **2,470** |                  0 |                        1 |                       2 |           0 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
..
.. The HiCUP categories same internal, re-ligation and contiguous corresponds to Diachromatic's un-ligated category.
.. HiCUP's same circularised category corresponds to the self-ligated category.
.. The wrong size category is corresponds to the sum of too short and too large chimeric fragments.
.. 99% of HiCUP's valid read pairs are also categorized as valid within Diachromatic.


Running Diachromatic's *align* subcommand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the alignment step: ::

    $ java -jar target/Diachromatic.jar align \
        -b /usr/bin/bowtie2 \
        -i /data/bt_indices/hg38 \
        -q prefix.truncated_R1.fq.gz \
        -r prefix.truncated_R2.fq.gz \
        -d hg38_DpnII_DigestedGenome.txt


The table lists all possible arguments:

+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| Short option | Long option                  | Example                                     | Required | Description                                                          | Default |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -q           | \\-\\-fastq-r1               | prefix.truncated_R1.fq.gz                   | yes      | Path to the truncated forward FASTQ file.                            | --      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -r           | \\-\\-fastq-r2               | prefix.truncated_R2.fq.gz                   | yes      | Path to the truncated forward FASTQ file.                            | --      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -b           | \\-\\-bowtie2                | /tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2 | yes      | Path to bowtie2 executable.                                          | --      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -i           | \\-\\-bowtie2-index          | /data/indices/bowtie2/hg38/hg38             | yes      | Path to bowtie2 index of the corresponding genome.                   | --      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -d           | \\-\\-digest-file            | /data/GOPHER/hg38_DpnII_DigestedGenome.txt  | yes      | Path to the digest file produced with GOPHER.                        | --      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -o           | \\-\\-out-directory          | cd4v2                                       | no       | Directory containing the output of the align subcommand.             | results |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -x           | \\-\\-out-prefix             | stim_rep1                                   | no       | Prefix for all generated files in output directory.                  | prefix  |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -p           | \\-\\-thread-num             | 15                                          | no       | Number of threads used by bowtie2.                                   | 1       |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -j           | \\-\\-output-rejected        | --                                          | no       | If set, a BAM file containing the reject read pairs will be created. | false   |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -l           | \\-\\-lower-frag-size-limit  | 50                                          | no       | Lower threshold for the size of sheared fragments.                   | 50      |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -u           | \\-\\-upper-frag-size-limit  | 1000                                        | no       | Upper threshold for the size of sheared fragments.                   | 1000    |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+
| -s           | \\-\\-self-ligtion-threshold | 3000                                        | no       | Upper threshold for the size of self-ligating fragments.             | 3000    |
+--------------+------------------------------+---------------------------------------------+----------+----------------------------------------------------------------------+---------+


Output files
~~~~~~~~~~~~

The default name of the BAM file containing all unique valid pairs that can be used for downstream analysis is:

    * ``prefix.valid_pairs.aligned.bam``


If ``--output-rejected`` is set, Diachromatic will output a second BAM file cointaing all rejected pairs:

    * ``prefix.rejected_pairs.aligned.bam``


Diachromatic uses optional fields of the SAM records to indicate the read pair category:

    * Un-ligated due to size (Tag: ``UL``)
    * Un-ligated due to same digest (Tag: ``ULSI``)
    * Self-ligated due to size (Tag: ``SL``)
    * Self-ligated due to same digest (Tag: ``SLSI``)
    * Too short chimeric  (Tag: ``TS``)
    * Too long chimeric  (Tag: ``TL``)
    * Valid pair (Tag: ``VP``)


Furthermore, there is an ``RO`` attribute that indicates the relative orientation of the pair:

    * Same strand forward: ``F1F2``, ``F2F1``
    * Same strand reverse: ``R1R2``, ``R2R1``
    * Inwards: ``F1R2``, ``F2R1``
    * Outwards: ``R2F1``, ``R1F2``


In addition, a file ``prefix.align.stats.txt`` is produced that contains summary statistics about the alignment step.


Finally, an R script ``prefix.frag.sizes.counts.script.R`` is generated that contains fragment size counts and can be
used to generate a plot as shown above.
In order to produce a PDF file, execute the script as follows: ::

    $ Rscript prefix.frag.sizes.counts.script.R

Or source the script from the R environment: ::


    > source("prefix.frag.sizes.counts.script.R")

