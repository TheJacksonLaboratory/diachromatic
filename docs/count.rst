Counting of valid read pairs between pairs of restriction fragments
===================================================================

Requirements regarding the input BAM file containing the valid pairs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you used Diachromatic with the subcommand ``align``, nothing has to be done. If the BAM file was produced in a
different way, make sure that the two reads of given pairs occur consecutively. Furthermore, make sure that duplicates
were previously removed.

GOPHER digest file.

BED file for active fragments.


Simple and twisted interaction counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

The interactions are written to a tab separated text file the followning name by default:

    * ``interaction.counts.table.tsv``

The structure of this file is similar to that of BED files. Each line stands for one pair of interacting fragments.
Consider the following example:

 ::

    chr7    42304777        42314850        A       chr3    152941166       152943990       I       12:2


The first three columns contain the coordinates of a restriction fragment on chromosome 7.
The ``A`` in column 4 indicates that this fragment is defined as active,
i.e. it is part of a viewpoint that was enriched using capture technology.
The information about active states of fragments originates either from the GOPHER digest file passed to Diachromatic
using the ``-d`` option or from the additional input file passed using the ``-a`` option.







