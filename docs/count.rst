Counting of valid read pairs between pairs of restriction fragments
===================================================================

Mapped Hi-C read pairs are typically transformed into contact matrices, whereby the pairs are counted between windows of
fixed size, typically 5 kbp (`Forcato et al., 2017 <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5493985/>`_) provide a review of
the methodology). Diachromatic was developed for *capture Hi-C*, which achieves
a much higher resolution than Hi-C. Therefore, for Diachromatic the read counts are determined for each
restriction digest.


Required input files
~~~~~~~~~~~~~~~~~~~~

GOPHER digest file
------------------

Due to the fact that the counts are determined on the restriction fragment level, the :ref:`digest file <rstdigest>` needs to be passed to ``Diachromatic count``. If the captured viewpoints were designed with GOPHER,
this file also includes information about active and inactive restriction fragments.




BAM file with unique valid pairs
--------------------------------

The second required input file contains the unique valid mapped read pairs in BAM format. If this file was generated using
Diachromatic with the :ref:`align <rstalign>` subcommand, nothing has to be done or taken care of. If the BAM file was produced in a different way,
make sure that the two reads of any given pair occur consecutively. Furthermore, make sure that duplicates were previously
removed.


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

    $ java -jar Diachromatic.jar count \
        -v prefix.valid_pairs.aligned.bam \
        -d hg38_DpnII_DigestedGenome.txt


+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| Short option | Long option          | Example                                                | Required | Description                                                      | Default |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -v           | --valid-pairs-bam    | prefix.valid_pairs.aligned.bam                         | yes      | Path to BAM file containing unique valid pairs.                  |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -d           | --digest-file        | /data/GOPHER/hg38_DpnII_DigestedGenome.txt             | yes      | Path to the digest file produced with GOPHER.                    |    --   |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -od          | --out-directory      | cd4v2                                                  | no       | Directory containing the output of the align subcommand.         | results |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+
| -op          | --out-prefix         | stim_rep1                                              | no       | Prefix for all generated files in output directory.              | prefix  |
+--------------+----------------------+--------------------------------------------------------+----------+------------------------------------------------------------------+---------+


Output files
~~~~~~~~~~~~

Interaction counts
------------------

The interactions are written to a tab separated text file that has the following name by default:

    * ``prefix.interaction.counts.table.tsv``

The structure of this file is similar to that of (`iBED <https://bioconductor.org/packages/release/bioc/vignettes/Chicago/inst/doc/Chicago.html#ibed-format-ends-with-ibed>`_) files. Each line stands for one pair of interacting fragments.
Consider the following example:

 ::

    chr7    42304777        42314850        A       chr7    152941166      152943990      I       12:2
    chr7    42304777        42314850        A       chr7    38624777       38625305       I       7:4

The first three columns contain the coordinates of a restriction fragment on chromosome 7. The ``A`` in column 4
indicates that this fragment is defined to be active, i.e. it is part of a viewpoint that was enriched using capture technology.
The information about active states of fragments originates from the GOPHER digest file passed to Diachromatic
using the ``-d`` option.

Read counts at interacting fragments
------------------------------------

Another file that is created contains the counts of reads at interacting fragments. By default the name of this file is:

    * ``prefix.interacting.fragments.counts.table.tsv``

The structure is again similar to that of BED files. Consider the following example:

 ::

    chr7    42304777       42314850        A       25
    chr7    152941166      152943990       I       14
    chr7    38624777       38625305        I       11

The first three columns contain the coordinates of interacting restriction fragments. This is again followed by either an ``A`` or ``I`` in column 4,
whereby ``A`` means active and ``I`` inactive. The fifth column contains the read counts aggregated from all
interactions that end in the corresponding fragment. For better understanding, compare these counts to the two
interactions given above.

Quality metrics
~~~~~~~~~~~~~~~

Fraction of singleton interactions (FSI)
----------------------------------------

It has been pointed out that the Cis/Trans ratio quality measure depends also on other factors such as the genome size and
number of chromosomes of the analyzed species (Wingett 2015). Diachromatic provides an alternative and possibly more robust quality metric that
can be used to access the extent of cross-ligation. Amongst the trans read pairs, we generally observe a large proportion
of restriction fragments that are connected by only a single read pair. The number of all possible different cross-ligation
events (including cis and trans) can roughly be estimated as the square number of all restriction fragments across the
entire genome. Given this huge number, we reasoned that it is very unlikely that the same cross-ligation event occurs
twice. Therefore, we defined the fraction of singleton interactions as the ratio of singleton read pairs and all read pairs.


Interaction count statistics
----------------------------

As for the other subcommands, a text file containing summary statistics is generated:

    * ``prefix.count.stats.txt``

This file contains:

    * The total number of processed read pairs.
    * The read pair counts broken down into the eight possible pair orientations.
    * Summary statistics about interactions between active and inactive fragments.
    * Quality metrics for experimental trouble shooting
        + Target Enrichment Coefficient (TEC): The fraction of reads that are mapped to active fragments.
        + Cross-ligation coefficient (CLC):	The fraction of trans read pairs.
        + Fraction of Singleton Interactions (FSI): The proportion of interactions consisting of only one read pair among all interactions.
            - This is an alternative quality metric that is intended to reflect the extend cross-ligation events.
