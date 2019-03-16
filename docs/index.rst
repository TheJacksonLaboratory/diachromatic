.. VPV documentation master file, created by
   sphinx-quickstart on Sun Sep 24 12:02:05 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Diachromatic's documentation!
========================================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   In silico genome digest <digest>
   Truncation of chimeric Hi-C reads <truncate>
   Mapping paired-end Hi-C reads <mapping>
   Counting unique valid pairs <count>
   Calling SNPs on Hi-C reads <allelspec>
   testing




Differential Analysis of Chromatin Interactions by Capture (Diachromatic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diachromatic_ is a Java application that implements a capture Hi-C preprocessing pipeline followed by analysis of differential chromatin interactions
("loopings").  Diachromatic is designed to work
with the capture probes as designed by `GOPHER <https://github.com/TheJacksonLaboratory/Gopher>`_ (see next section).

.. _Diachromatic: https://github.com/TheJacksonLaboratory/diachromatic


Generator Of Probes for capture Hi-C Expriments at high Resolution (GOPHER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GOPHER_ is a Java application designed to help design capture probes
for capture Hi-C and related protocols. Capture Hi-C (CHC) is based
on the Hi-C protocol but uses capture baits (similar to whole-exome sequencing) to enrich a set of viewpoints.
Commonly, the viewpoints represent proximal promoter regions (surrounding the transcription start site [TSS]) of
genes of interest or of all protein-coding genes. CHC detects interactions between viewpoint regions and distal
enhancers (or other genomic regions) and has been most commonly performed with the 4-cutter DpnII or with the 6-cutter
HindIII.

.. _GOPHER: https://github.com/TheJacksonLaboratory/Gopher


Quick start
~~~~~~~~~~~
Diachromatic requires Java 8 or higher to run. The source code of Diachromatic can be downloaded
from the Diachromatic `GitHub page <https://github.com/TheJacksonLaboratory/diachromatic>`_.

To build the application, clone the repository and create the Java app with maven. ::

    $ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
    $ cd diachromatic
    $ mvn package

To test whether the build process was successful, enter the following command: ::

    $ java -jar target/Diachromatic.jar

You should see a help message in the shell.


The mapping step of the diachromatic pipeline relies on
`bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_. If needed, install bowtie2
on your system. For instance, on Debian linux systems bowtie2 can be installed with the
following command. ::

  $ sudo apt-get install bowtie2


Preparation of the bowtie2 index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The prebuilt ``bowtie2`` indices for human hg19 (3.5 GB) and other genome builds can be downloaded from the
`bowtie2 website`_. After downloading the correct archived file to your computer, unpack it with: ::

    $ unzip hg19.zip
        Archive:  hg19.zip
        inflating: hg19.1.bt2
        inflating: hg19.2.bt2
        inflating: hg19.3.bt2
        inflating: hg19.4.bt2
        inflating: hg19.rev.1.bt2
        inflating: hg19.rev.2.bt2
        inflating: make_hg19.sh

In the following pages, we will call the path to the directory where
the index was unpacked **/path/to/bowtie2index/**. Substitute this with the actual path on your computer.

.. _bowtie2 website: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml


The remaining pages of this site explain the various steps of preprocessing and analyzing
capture Hi-C data with diachromatic.
