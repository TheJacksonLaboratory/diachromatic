
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

The first line of this file contains the column names, and all other lines correspomd to one restriction fragment.
Each line consists of 14 fields that are described in the table below. The features of the fields 8 to 14 will be used for the Poisson-regression based normalization scheme.

+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| Column | Name                    | Example  | Description                                                                                                                              |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 1      | Chromosome              | chr1     | Name of the reference sequence.                                                                                                          |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 2      | Fragment_Start_Position | 18376    | 1-based start position of the restriction fragment.                                                                                      |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 3      | Fragment_End_Position   | 18392    | 1-based end position of the restriction fragment.                                                                                        |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 4      | Fragment_Number         | 42       | Consecutive fragment number.                                                                                                             |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 5      | 5'_Restriction_Site     | DpnII    | Name of the enzyme responsible for the cut at the 5' end of the fragment.                                                                |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 6      | 3'_Restriction_Site     | DpnII    | Name of the enzyme responsible for the cut at the 3' end of the fragment. May be different from field 5 if more than one enzyme is used. |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 7      | Length                  | 1245     | Length of the fragment.                                                                                                                  |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 8      | 5'_GC_Content           | 0.500    | GC content of the upstream margin (GOPHER's default margin size is 250 bp).                                                              |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 9      | 3'_GC_Content           | 0.500    | GC content of the downstream margin.                                                                                                     |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 10     | 5'_Repeat_Content       | 0.138    | Repeat content of the upstream margin.                                                                                                   |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 11     | 3'_Repeat_Content       | 0.126    | Repeat content of the downstream margin.                                                                                                 |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 12     | Activation status       | T (or F) | Activation status of the fragment. A fragment is active, if it has at least one probe.                                                   |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 13     | 5'_Probes               | 2        | Number of probes for the upstream margin.                                                                                                |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
| 14     | 3'_Repeat_Content       | 0        | Number of probes for the downstream margin.                                                                                              |
+--------+-------------------------+----------+------------------------------------------------------------------------------------------------------------------------------------------+
