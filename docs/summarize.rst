
Summarize command
=================


.. Quality metrics
.. ~~~~~~~~~~~~~~~
.. Diachromatic outputs a text file with the quality metrics for each run. The following texts provides possible interpretations
.. of these quality metrics and exemplary numbers for the `CTCF depletion Hi-C datasets of Nora et al. 2017`_.
..
.. .. _CTCF depletion Hi-C datasets of Nora et al. 2017: https://www.ncbi.nlm.nih.gov/pubmed/2852575
..
..
.. Percentage of uniquely mapped pairs
.. -----------------------------------
..
.. Percentage of truncated input read pairs that were paired, i.e. both reads were uniquely mapped to the genome.
.. For the CTCF depletion datasets, percentages range from 48.54% to 56.66%.
..
..
.. Hi-C pair duplication rate (HPDR)
.. ---------------------------------
..
.. For Hi-C, the removal of duplicates must take into account the chimeric nature of the underlying fragments.
.. The HPDR is defined as the percentage of uniquely mapped pairs that were removed because they were recognized to be *Hi-C duplicates*.
.. Usually, high duplication rates indicate sequencing libraries with low complexity.
.. For the CTCF depletion data, the proportion of unique read pairs amongst all uniquely mapped read pairs ranges between
.. 1.26% and 21.13%.
..
..
.. Percentages of different read pair categories
.. ---------------------------------------------
..
.. The categorization scheme subdivides the set of all uniquely mapped unique pairs into disjoint subsets.
.. The percentages of the individual categories may be useful for experimental troubleshooting.
..
.. **Percentage of un-ligated read pairs:** High percentages of un-ligated pairs indicate poor enrichment for ligation junctions, i.e. the streptavidin pull-down of biotinylated Hi-C ligation did not perform well. For the CTCF depletion data, the percentages of un-ligated pairs range between 9.61% and 26.17%.
..
.. **Percentage of self-ligated read pairs:** In practice, self-ligation seems not to occur  very often. For the CTCF depletion data, the percentages of self-ligated pairs range between 0.95% and 1.89%.
..
.. **Percentage of too short chimeric read pairs:** A high percentage (5%<) of too short chimeric fragments may indicate that the chosen lower size threshold for sheared fragments (``-l``) does not match the experimental settings. Diachromatic generates a plot for distribution of fragment sizes (see below) may provide guidance.
..
.. **Percentage of too large chimeric read pairs:** Essentially, the same applies as for the too short chimeric category.
..
.. **Percentage of valid read pairs:** The more, the better. For the the CTCF depletion data, percentages range between 62.30% and 81.35%.
..
..
.. Yield of valid pairs (YVP)
.. --------------------------
..
.. Percentage of truncated input read pairs that  are not
.. categorized as artifactual by any of the quality control steps, and therefore can be used for downstream analysis.
.. The YVP reflects the overall efficiency of the Hi-C protocol.
.. For the the CTCF depletion data, the percentages range between 24.37% and 42.77%.
..
..
.. Cross-ligation coefficient (CLC)
.. --------------------------------
..
.. Valid read pairs arising from genuine chromatin-chromatin interactions between different chromosomes cannot be
.. distinguished from those arising from **cross-ligation** events.
.. Based on the assumption that random cross-ligations between DNA fragments of different chromosomes (*trans*) occur more
.. likely as compared to cross-ligations between DNA fragments of the same chromosome (*cis*), the ratio of the numbers of cis
.. and trans read pairs is taken as an indicator of poor Hi-C libraries (Wingett 2015, Nagano 2015).
.. Within Diachromatic, the CLC is calculated as the proportion of unique valid trans read pairs amongst all unique valid read pairs.
.. For the CTCF depletion dataset, percentages range between 18.48% and 28.24%.
..
..
.. Re-ligation coefficient (RLC)
.. -----------------------------
..
.. Percentage of uniquely mapped unique pairs that did not arise from fragments with dangling-ends, i.e. ends that correspond
.. to un-ligated restriction enzyme cutting sites.
.. The RLC is intended to reflect the efficiency of the re-ligation step
.. and could possibly be used to improve experimental
.. conditions for re-ligation.
.. For the CTCF depletion dataset, percentages range between 97.04% and 98.92%.
..
..
.. Size distribution of chimeric and un-ligated fragments
.. ------------------------------------------------------
..
.. The plot of fragment size distributions is intended to serve as a kind of sanity check.
.. Deviations from bell-shaped curve progressions should be thoroughly scrutinized.
.. Furthermore, the plot might be useful for the adjustment of Diachromatic's size thresholds T1\ :sub:`min` and T1\ :sub:`max`.
.. For instance, a high number of read pairs that are categorized as *too large* could indicate that the actual size of
.. sheared fragments is larger on average.
.. In such cases, the plot can be used to choose good thresholds.
..
.. For the size distribution of chimeric fragments (**black**), the chimeric sizes of all read pairs that were categorized
.. as either as *valid*, *too short* or *too long* are determined.
.. Enriched chimeric fragments (**red**) form a subset of all chimeric fragments, whereby either the read R1 or R2 is assigned
.. to a digest that is flagged as selected in the digest file passed to Diachromatic.
.. For the size distribution of un-ligated fragments (**blue**) the distances between all inward pointing read pairs mapping
.. to the same chromosome (*cis*) are determined.
..
.. .. figure:: img/size_distribution_plot.png
..     :align: center

