Counting of valid read pairs between pairs of restriction fragments
===================================================================

Mapped Hi-C read pairs are typically transformed into contact matrices, whereby the pairs are counted between windows of
fixed size, typically 5 kbp (ref). Diachromatic was developed in the first place for *capture Hi-C* data which achieves
a much higher resolution as compared to Hi-C. Therefore, for Diachromatic the read counts are determined on the
restriction fragment level.


Required input files
~~~~~~~~~~~~~~~~~~~~

GOPHER digest file
------------------

Due to the fact that the counts are determined on the restriction fragment level, the digest file `initially produced
using GOPHER`_ needs to be passed to ``Diachromatic count``. If the captured viewpoints were designed with GOPHER,
this file also includes information about active and inactive restriction fragments. If this is not the case,
you can pass an additional BED file containing the coordinates of active digests. Note that the coordinates must exactly
match the coordinates in the digest file. The coordinates in GOPHER's digest file are 1-based. Therefore, the BED file
for active fragments does strictly speaking not comply with the BED format that has 0-based start coordinates and
1-based end coordinates.

.. _initially produced using GOPHER: digest.html


BAM file with unique valid pairs
--------------------------------

The second required input file contains the unique valid mapped read pairs in BAM format. If this file was generated `using
Diachromatic with the align subcommand`_, nothing has to be done or taken care of. If the BAM file was produced in a different way,
make sure that the two reads of any given pair occur consecutively. Furthermore, make sure that duplicates were previously
removed.

.. _using Diachromatic with the align subcommand: mapping.html


Simple and twisted read pairs and counts of directed interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diachromatic aggregates read pairs whose 5' end positions map to the same pair of restriction fragments into interaction counts,
whereby one special feature of Diachromatic is that the relative orientation of read pairs is taken into account.
Inward and outward pointing read pairs (F1R2,F2R1,R1F2,R2F1) are referred to as **simple**, whereas pairs for which the two reads are
pointing in the same direction (R1R2,R2R1,F1F2,F2F1) are referred to as **twisted**.
Therefore, a given pair of interacting restriction fragments is assigned two interaction counts separated by a colon
character. For instance, if we have ``12:2`` for a given pair of restriction fragments, this means that there are ``12``
simple and ``2`` read pairs. The two counts can be used in order to distinguish **directed** from **undirected** interactions
(see manuscript). For the given example ``12:2``, could be considered as directed.
However, if the same number of read pairs were distributed like ``6:8``, the interaction could be considered as undirected.
At the moment, Diachromatic does not provide any rules or statistical framework in order to distinguish directed from
undirected interactions.


Running Diachromatic's truncation subcommand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the counting step. ::

    $ java -jar Diachromatic.jar count -v prefix.valid_pairs.aligned.bam -d hg38_DpnII_DigestedGenome.txt


+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| Short option | Long option          | Example                                                | Required | Description                                                      | Default |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -v           | --valid-pairs-bam    | prefix.valid_pairs.aligned.bam                         | yes      | Path to BAM file containing unique valid pairs.                  |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -d           | --digest-file        | /data/GOPHER/hg38_DpnII_DigestedGenome.txt             | yes      | Path to the digest file produced with GOPHER.                    |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -a           | --active-digest-file | /data/GOPHER/hg38_DpnII_active_digests_cd4v2_genes.bed | no       | Path to a BED file containing the coordinates of active digests. |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -od          | --out-directory      | cd4v2                                                  | no       | Directory containing the output of the align subcommand.         | results |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -op          | --out-prefix         | stim_rep1                                              | no       | Prefix for all generated files in output directory.              | prefix  |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+


Output files
~~~~~~~~~~~~

The interactions are written to a tab separated text file that has the following name by default:

    * ``interaction.counts.table.tsv``

The structure of this file is similar to that of iBED files. Each line stands for one pair of interacting fragments.
Consider the following example:

 ::

    chr7    42304777        42314850        A       chr3    152941166       152943990       I       12:2

The first three columns contain the coordinates of a restriction fragment on chromosome 7. The ``A`` in column 4
indicates that this fragment is defined to be active, i.e. it is part of a viewpoint that was enriched using capture technology.
The information about active states of fragments originates either from the GOPHER digest file passed to Diachromatic
using the ``-d`` option or from the additional input file passed using the ``-a`` option.
