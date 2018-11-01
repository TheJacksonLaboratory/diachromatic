
Creating an in silico restriction digest map for Diachromatic using GOPHER
==========================================================================

The Capture Hi-C (CHC) protocol involves the restriction digestion of a sample and the downstream analysis assigns
reads to pairs of restriction fragments. Therefore, the ``align`` subcommand of Diachromatic requires a list of all
restriction fragments that result from the in silico digestions of a given genome with the chosen enzyme or enzymes.
Such lists can be generated using the GOPHER_ software for design of capture Hi-C viewpoints and probes. The TSV
formatted file exported from GOPHER can be passed to Diachromatic using the ``-m`` or ``digest`` option.

.. _GOPHER: https://github.com/TheJacksonLaboratory/Gopher

Performing the in silico digest using GOPHER
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the GOPHER software was used for the design of the given capture Hi-C experiment, you can just open the corresponding
GOPHER project and export the required file via the export menu. In such cases the exported file also includes the
information about active and inactive digests, i.e. if there are probes associated with given digests.

.. figure:: img/output_export.png

If you did not perform the design using GOPHER, you will have to `setup a new project`_.
It's sufficient to specify the parameters within the *Data sources* section and to selected used enzyme within the
*Design parameters* section. If you prepare the digest map in this way, all digests will be marked as inactive.

.. _setup a new project: https://gopher.readthedocs.io/en/latest/02_gui_data.html

.. figure:: img/digest_parameters.png

In order to retrospectively mark digests as active an additional BED file can be passed to Diachromatic using the
option ``-a``. This file should contain all digests that are desired to be marked as active. The coordinates of the digests
have to be exactly the same as in the digest file produced by GOPHER. It is useful to use ``bedtools intersect`` in order
to make sure that the coordinates of the digests are identical between the two files.

Format of the digest file
~~~~~~~~~~~~~~~~~~~~~~~~~


+--------+-------------------------+----------+
| Column | Name                    | Example  |
+--------+-------------------------+----------+
| 1      | Chromosome              | chr1     |
+--------+-------------------------+----------+
| 2      | Fragment_Start_Position | 18376    |
+--------+-------------------------+----------+
| 3      | Fragment_End_Position   | 18392    |
+--------+-------------------------+----------+
| 4      | Fragment_Number         | 42       |
+--------+-------------------------+----------+
| 5      | 5'_Restriction_Site     | DpnII    |
+--------+-------------------------+----------+
| 6      | 3'_Restriction_Site     | DpnII    |
+--------+-------------------------+----------+
| 7      | Length                  | 1245     |
+--------+-------------------------+----------+
| 8      | 5'_GC_Content           | 0.500    |
+--------+-------------------------+----------+
| 9      | 3'_GC_Content           | 0.500    |
+--------+-------------------------+----------+
| 10     | 5'_Repeat_Content       | 0.138    |
+--------+-------------------------+----------+
| 11     | 3'_Repeat_Content       | 0.126    |
+--------+-------------------------+----------+
| 12     | Selected                | F (or T) |
+--------+-------------------------+----------+
| 13     | 5'_Probes               | 2        |
+--------+-------------------------+----------+
| 14     | 3'_Repeat_Content       | 0        |
+--------+-------------------------+----------+


+----------------+----------+---------+--------+----------+----------+
| Chromosome     |Start_Pos | End_Pos | Frag_N | 5'_RSite | 3'_RSite |
+================+==========+=========+========+==========+==========+
| chrUn_gl000243 | 1        |  1816   | 1      | None     | HindIII  |
+----------------+----------+---------+--------+----------+----------+
| chrUn_gl000243 | 1817     |   4040  | 2      | HindIII  | HindIII  |
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




Data for Poisson regression scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We are extending the table to include data that we will use for the Poisson-regression based normalization scheme.



