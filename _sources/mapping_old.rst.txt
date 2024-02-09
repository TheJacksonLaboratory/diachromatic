
Categorization of read pairs
----------------------------

.. **1. Un-ligated:** The read pair points inwards and the distance between the two 5' end positions d\ :sub:`u` is smaller than T1\ :sub:`max` or both reads map to the same digest.
..
.. **2. Self-ligated:** The read pair points outwards and the calculated size of self-ligating fragments d\ :sub:`s` is smaller than a predefined self-ligation threshold T2\ :sub:`max` (Default: 3000) or both reads map to the same digest.
..
.. **3. Short chimeric:** The read pair is not in the un-ligated or self-ligated category and the calculated size d\ :sub:`c` is smaller than a specified lower threshold threshold T1\ :sub:`min` (Default: 50).
..
.. **4. Long chimeric:** The read pair  is not in the un-ligated or self-ligated category and the calculated size d\ :sub:`c` is greater than a specified lower threshold T1\ :sub:`max` (Default: 800).
..
.. **5. Valid (chimeric):** All remaining chimeric read pairs.
..
..
.. The decision as to whether a read-pair is valid or not is made according to:
..
.. **1.** Read pairs that map to different chromosomes or to the same strand cannot originate from un-ligated or self-ligated fragments. Therefore, they are categorized as chimeric read pairs that are valid, if the size d\ :sub:`s` is within the specified range.
..
.. **2.** Read pairs that point inwards might originate from un-ligated fragments. In such cases, the distance between the 5' end positions of the mapped reads d\ :sub:`u` corresponds to the size of the  sequenced fragment. In order to assign read pairs to the un-ligated category, we use an upper size threshold T\ :sub:`1` that should reflect the maximum plausible size of sheared fragments. Furthermore, inward pointing read pairs that map to the same digest are categorized as un-ligated.
..
.. **3.** Read pairs that point outwards might originate from self-ligated fragments. In such cases, the size d\ :sub:`s` of the potentially underlying self-ligated fragment is calculated as described above, and compared to an upper size threshold T\ :sub:`2` for self-ligated fragments. Outward pointing read pairs with d\ :sub:`s` smaller than T\ :sub:`2` are assigned to the self-ligated category. Furthermore, outward pointing read pairs that map to the same digest are categorized as self-ligated.
..
.. **4.** Read pairs arising from chimeric fragments (not un- or self-ligated) are further distinguished. Read pairs with size d\ :sub:`s` outside the specified size range of sheared fragments will be categorizesd as too small or too large, and all remaining read pairs are categorized as valid.



.. Dangling end read pairs
.. -----------------------
..
.. Fragment ends that corresponding to restriction enzyme cutting sites are referred to as dangling ends.
.. In theory, fragments of all categories may have dangling ends. Therefore, there is no separate class for dangling ends.
.. However, the number of dangling end read pairs within each of the five disjoint categories is determined and reported.



.. Dichromatic vs. HiCUP categories
.. --------------------------------
..
.. When HiCUP is executed with the ``--keep`` flag, it will create a directory containing BAM files for the individual read pair
.. categories. We applied HiCUP to the associated test data, converted the BAM files back to FASTQ format
.. and applied Diachromatic to the FASTQ files.
..
.. The following table shows the numbers of read pairs within the categories of HiCUP and Diachromatic.
..
.. For instance, HiCUP categorized 13,760 read pairs as *same internal* and 13,722 of these are uniquely mapped using Diachromatic.
.. The small differences between these numbers may be due to different bowtie versions or settings.
.. 13,645 of these uniquely mapped read pairs are categorized as un-ligated, which is the correct category for those read pairs
.. because according to our logic *same internal* read pairs correspond to un-ligated fragments.
.. However, in total 77 *same internal* read pairs are categorized as *chimeric* read pairs, which is contradictory.
.. Further investigation revealed that the 5' end positions of those read pairs are indeed mapped to the same digest but also to the same strand.
.. According to the logic implemented in Diachromatic (see decision tree) read pairs mapped to the same strand are automatically categorized as chimeric,
.. because the concept of the Hi-C fragment formation cannot explain such read pairs.
..
.. The next HiCUP category is *re-ligation*. For Diachromatic, all 1060 read pairs are mapped uniquely.
.. 58 read pairs are categorized as *self-ligated*. Further investigation of these read pairs revealed that all pairs
.. are outward pointing, which is correct for *self-ligated* pairs.
.. The 5 *re-ligation* read pairs that are categorized as *chimeric too short* are outward pointing as well but d\ :sub:`u` is greater than
.. the self-ligation threshold. However, the calcluated size calculated d\ :sub:`c` is smaller than lower threshold for sheared fragments.
..
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **HiCUP** (rows) vs. **Diachromatic** (columns)      | **# Processed pairs** | **# Uniquely mapped pairs** | **# Un-ligated** | **# Self-ligated** | **# Chimeric too short** | **# Chimeric too long** | **# Valid** |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same internal**                                    |                13,760 |                      13,722 |       **13,645** |                  0 |                       13 |                      39 |          25 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Re-ligation**                                      |                 1,060 |                       1,060 |          **842** |                 58 |                        5 |                      49 |         106 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Contiguous**                                       |                    58 |                          58 |           **53** |                  0 |                        1 |                       0 |           4 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same circularised**                                |                   428 |                         428 |                3 |            **425** |                        0 |                       0 |           0 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Wrong size**                                       |                10,321 |                      10,267 |                2 |                  0 |                **1,003** |               **9,181** |          81 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Valid**                                            |                25,915 |                      25,851 |                1 |                  5 |                      290 |                       6 |  **25,549** |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. |                                                      |                       |                             |                  |                    |                          |                         |             |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
.. | **Same dangling ends**                               |                 2,475 |                       2,473 |        **2,470** |                  0 |                        1 |                       2 |           0 |
.. +------------------------------------------------------+-----------------------+-----------------------------+------------------+--------------------+--------------------------+-------------------------+-------------+
..
.. The HiCUP categories same internal, re-ligation and contiguous corresponds to Diachromatic's un-ligated category.
.. HiCUP's same circularised category corresponds to the self-ligated category.
.. The wrong size category is corresponds to the sum of too short and too large chimeric fragments.
.. 99% of HiCUP's valid read pairs are also categorized as valid within Diachromatic.

