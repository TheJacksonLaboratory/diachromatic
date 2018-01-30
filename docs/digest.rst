
Diachromatic: Creating an in silico restriction digest map of the genome
========================================================================

The Capture Hi-C (CHC) protocol involves the restriction digestion of a sample
and the downstream analysis assigns reads to pairs of restriction fragments
(See todo for details on the CHC protocol). Therefore, the first step of
the diachromatic analysis pipeline is the preparation of an in silico digest
of the genome corresponding to the sample and the restriction enzyme used in
the experiment.

Performing the digest
~~~~~~~~~~~~~~~~~~~~~
The command to perform in silico digestion is: ::

    $ java -jar Diachromatic.jar digest -g <path> -e <enzyme> [-o <outfile>]

The meaning of the options is
   * -g <path> The path to the directory with genome FASTA files (one per chromosome; for instance, the genome fasta files downloaded from the UCSC Genome Browser are in this format). The FASTA files need to be unzipped.
   * -e <enzyme> The symbol of the restriction enzyme used in the Capture Hi-C experiment.
   * [-o <outfile>] Name and path of the output file. This flag is optional and if it is not passed, the default name of ``diachromaticDigest.txt`` will be used.

For example, the following command will digest the genome files found in the directory ``hg19`` using
the restriction enzyme ``HindIII`` and will produce an output file called ``hg19HindIIIdigest.txt``. ::


   $ java -jar Diachromatic.jar digest -g /path/to/hg19/ -e HindIII -o hg19HindIIIdigest.txt

Output file format
~~~~~~~~~~~~~~~~~~
The output file has the following format.


+----------------+----------+---------+--------+----------+----------+
| Chromosome     |Start_Pos | End_Pos | Frag_N | 5'_RSite | 3'_RSite |
+================+==========+=========+========+==========+==========+
| chrUn_gl000243 | 1        |  1816   | 1      | None     | HindIII  |
+----------------+----------+---------+--------+----------+----------+
| chrUn_gl000243 | 1817     |   4040  |   2    | HindIII  | HindIII  |
+----------------+----------+---------+--------+----------+----------+
| chrUn_gl000243 | 4041     | 10052   | 3      | HindIII  | HindIII  |
+----------------+----------+---------+--------+----------+----------+
| chrUn_gl000243 | 10053    |   10364 | 4      | HindIII  | HindIII  |
+----------------+----------+---------+--------+----------+----------+
| chrUn_gl000243 | 10401    |   10698 | 5      | HindIII  | HindIII  |
+----------------+----------+---------+--------+----------+----------+

The file shows all restriction fragments with their chromosomes, start and end positions,
and the Restriction enzyme responsible for the 5' and 3' cut. This file will be used in
the ``map`` step of the Diachromatic pipeline.




**ToDo** offer option to restrict to canonical chromosomes.

**To Do** Check whether the first entry is correct for DpnII (not shown now). If I run the truncate command with DpnII I get

+-------------+-----------+---------+--------+----------+----------+
| Chromosome  | Start_Pos | End_Pos | Frag_N | 5'_RSite | 3'_RSite |
+=============+===========+=========+========+==========+==========+
| chr22       | 1         | 16050000| 1      | None     | DpnII    |
+-------------+-----------+---------+--------+----------+----------+
| chr22       | 16050001  | 16050155| 2      | DpnII    | DpnII    |
+-------------+-----------+---------+--------+----------+----------+
| chr22       | 16050156  | 16050444| 3      | DpnII    | DpnII    |
+-------------+-----------+---------+--------+----------+----------+

and the first entry is correct.



