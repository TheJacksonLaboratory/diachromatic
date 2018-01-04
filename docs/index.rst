.. VPV documentation master file, created by
   sphinx-quickstart on Sun Sep 24 12:02:05 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Diachromatic's documentation!
===============================
.. toctree::
	:maxdepth: 2
	:caption: Contents:

	CHC

Diachromatic
~~~~~~~~~~~~
Diachromatic (**Di**fferential **A**nalysis of **Chromat**in **I**nteractions by **C**apture)
implements a capture Hi-C preprocessing pipeline followed by analysis of differential chromatin interactions
("loopings").  Diachromatic is a Java application. The preprocessing pipeline is based on the
Perl scripts of HiCUP and produces comparable results. Diachromatic is designed to work
with the capture probes as designed by VPV (todo--link to RTD site).



ViewPointViewer (VPV)
~~~~~~~~~~~~~~~~~~~~~~~~~

VPV is a Java application designed to help design capture probes
for capture Hi-C and related protocols. Capture Hi-C (CHC) is based
on the Hi-C protocol but uses capture baits (similar to whole-exome sequencing) to enrich a set of viewpoints.
Commonly, the viewpoints represent proximal promoter regions (surrounding the transcription start site [TSS]) of
genes of interest or of all protein-coding genes.

- CHC detects interactions between viewpoint regions and distal enhancers (or other genomic regions).
- CHC has been most commonly performed with the 4-cutter DpnII or with the 6-cutter HindIII.
- For more information, see the VPV GitHub page at https://github.com/TheJacksonLaboratory/VPV.



Quick start
~~~~~~~~~~~~~~~~~~~~~~~~~
Diachromatic requires Java 8 or higher to run. The source code of Diachromatic can be downloaded
from the Diachromatic GitHub repository and the application can be built using maven (see the GitHub page for instructions).

This site provides detailed explanations and tips for the various steps of preprocessing and analyzing
capture Hi-C data with diachromatic.
 



