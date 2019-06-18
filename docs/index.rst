
Welcome to Diachromatic's documentation
=======================================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Program setup <setup>
   Tutorial with test data <tutorial>
   In silico genome digest <digest>
   Truncation of chimeric Hi-C reads <truncate>
   Mapping paired-end Hi-C reads <mapping>
   Counting unique valid pairs <count>
   Summarize results <summarize>

..   Calling SNPs on Hi-C reads <allelspec>


Differential Analysis of Chromatin Interactions by Capture (Diachromatic)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Diachromatic_ is a Java application that implements a capture Hi-C preprocessing pipeline followed by analysis of differential chromatin interactions ("loopings"). Diachromatic is designed to work with the capture probes as designed by `GOPHER <https://github.com/TheJacksonLaboratory/Gopher>`_.

.. _Diachromatic: https://github.com/TheJacksonLaboratory/diachromatic


Generator Of Probes for capture Hi-C Expriments at high Resolution (GOPHER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GOPHER_ is a Java application designed to help design capture probes for capture Hi-C and related protocols. Capture Hi-C (CHC) is based on the Hi-C protocol but uses capture baits (similar to whole-exome sequencing) to enrich a set of viewpoints. Commonly, the viewpoints represent proximal promoter regions (surrounding the transcription start site [TSS]) of genes of interest or of all protein-coding genes. CHC detects interactions between viewpoint regions and distal enhancers (or other genomic regions) and has been most commonly performed with the 4-cutter DpnII or with the 6-cutter HindIII.

.. _GOPHER: https://github.com/TheJacksonLaboratory/Gopher

