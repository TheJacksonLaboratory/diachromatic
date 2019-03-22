
Truncation of chimeric reads
============================

Ligation junctions, chimeric fragments and chimeric reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Valid Hi-C read pairs originate from chimeric fragments with DNA from two different loci linked by the ligation junction.
Most experimental protocols for capture Hi-C fill in the sticky ends with biotinylated nucleotides.
This alters the original restriction enzyme motif in a characteristic way.
For instance, the enzyme HinDIII cuts the DNA in the following manner:

.. figure:: img/HindIII_cut_site.png
    :align: center

Therefore, if some sequence contains a HindIII site, '5-NAAGCTTN-3', restriction followed by filling in and religation of the resulting blunt ends with produce the following hybrid sequence '5-NAAGCTAGCTTN-3'.
Diachromatic searches for hybrid patterns like this according the restriction enzyme or enzymes that were used, and truncates the reads accordingly.
Rarely, protocols do not fill in the ends but instead perform religation and thereby recreate the original restriction motif.

.. figure:: img/sticky_and_blunt_ends.png
    :align: center

The sonication step of the Hi-C protocol may introduce breakpoints near restriction enzyme cutting
sites. If the breakpoint occurs at a distance smaller than one read length, this will result in a chimeric read that
cannot be mapped to the reference sequence.
The truncation step of the pipeline attempts to address this problem by deleting the sequence that is downstream of
the enzyme recognition site.

.. figure:: img/chimeric_reads.png
    :align: center



Running Diachromatic's *truncate* subcommand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the truncation step. ::

    $ java -jar Diachromatic.jar truncate \
        -q test1.fastq \
        -r test2.fastq \
        -e HindIII \
        -x prefix \
        -o outdir


The following table lists all possible arguments:

+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| Short option | Long option     | Example       | Required | Description                                              | Default |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -q           | --fastq-r1      | forward.fq.gz | yes      | Path to the forward FASTQ file.                          |    --   |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -r           | --fastq-r2      | reverse.fq.gz | yes      | Path to the reverse FASTQ file.                          |    --   |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -e           | --enzyme        | HindIII       | yes      | Symbol of the restriction enzyme.                        | null    |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -s           | --sticky-ends   | false         | no       | True, if no fill-in of sticky ends was performed.        | false   |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -o           | --out-directory | cd4v2         | yes      | Directory containing the output of the truncate command. | results |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -x           | ---out-prefix   | stim_rep1     | yes      | Prefix for all generated files in output directory.      | prefix  |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+

Output files
~~~~~~~~~~~~

The default names of the truncated and gzipped FASTQ files are:

    * ``prefix.truncated_R1.fastq.gz``
    * ``prefix.truncated_R2.fastq.gz``

In addition, a file is produced that contains summary statistics about the truncation step.

    * ``prefix.truncation.stats.txt``

