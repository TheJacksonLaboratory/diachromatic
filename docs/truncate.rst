
Truncation of chimeric reads
============================

Ligation junctions, chimeric fragments and chimeric reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Valid Hi-C read pairs stem from chimeric fragments consisting of DNA from two different loci linked by the ligation
junction. Depending on whether the sticky ends of the dangling ends were filled or not, the ligation junction consist
of either one or two restriction enzyme cutting motifs.

.. figure:: img/sticky_and_blunt_ends.png
    :align: center

The sonication step of the Hi-C protocol may introduce breakpoints near restriction enzyme cutting
sites. If the breakpoint occurs at a distance smaller than one read length, this will result in a chimeric read that
cannot be mapped to the reference sequence.

.. figure:: img/chimeric_reads.png
    :align: center

The truncation step of the pipeline attempts to address this situation by deleting the sequence that is downstream of
the enzyme recognition site.

Running Diachromatic's truncation subcommand
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the truncation step. ::

    $ java -jar Diachromatic.jar truncate -q test1.fastq -r test2.fastq -e HindIII


The following table lists all possible arguments.

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
| -od          | --out-directory | cd4v2         | no       | Directory containing the output of the truncate command. | results |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+
| -op          | ---out-prefix   | stim_rep1     | no       | Prefix for all generated files in output directory.      | prefix  |
+--------------+-----------------+---------------+----------+----------------------------------------------------------+---------+

Output files
~~~~~~~~~~~~

The default names of the truncated and gzipped FASTQ files are:

    * prefix.truncated_R1.fastq.gz
    * prefix.truncated_R2.fastq.gz

In addition, a file is produced that contains summary statistics about the truncation step.

    * prefix.truncation.stats.txt

