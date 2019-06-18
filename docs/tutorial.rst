
Tutorial
========

This tutorial shows how to use Diachromatic for processing and quality control of Capture Hi-C reads.


Test dataset
~~~~~~~~~~~~

To get the data, visit this [ftp server](ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/) and download the two test read files or use: ::

	wget ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/test*.fastq


Truncation
~~~~~~~~~~

The first step of processing with data with Diachromatic is processing of the raw FASTQ files to recognize and truncate
reads with filled-in ligation juctions, which indicate reads that include the junction of the chimeric CHC fragment. If
a single read is chimeric, it is not possible to map it to one locus, and therefore the 3' portion of the chimeric read
is removed ("truncated"), leaving behind the 5' portion of the read that should map to a specific locus. If the 5' sequence
is too short to be mapped, the entire read pair is discarded. This is performed with the truncate command. ::

    $ java -jar Diachromatic.jar truncate \
        -q test_1.fastq \
        -r test_2.fastq \
        -e HinDIII \
        -x prefix \
        -o outdir


In practice, only about XXXX percent of the readpairs are truncated.


Mapping
~~~~~~~

The second step of the pipeline is to map the truncated read pairs to the target genome. Diachromatic uses bowtie2 to perform the
mapping, and then creates a BAM file containing the valid read pairs. If desired, Diachromatic also outputs BAM files
with the discarded (arterfactual or unmappable reads).

To use the mapping file, we need to have a file that shows the locations of restriction digests across the genome.
Diachromatic users can use GOPHER to create probes and the digest file TODO -- MORE DETAIL.

It is recommended that users download bowtie index files from the bowtie site. In the following, we have
downloaded and extracted these files to a directory called ``/data/bt_indices`` (in this directory, there are multiple index files
hg19.1.bt2, hg19.3.bt2, hg19.rev.1.bt2, hg19.2.bt2, hg19.4.bt2, hg19.rev.2.bt2).

Use the following command to run the alignment step. ::

    $ java -jar target/Diachromatic.jar align \
        -b /usr/bin/bowtie2 \
        -i /data/bt_indices/hg19 \
        -q prefix.truncated_R1.fastq.gz \
        -r prefix.truncated_R2.fastq.gz \
        -d hg19_HinDIII_DigestedGenome.txt \
        -x prefix \
        -o outdir


Counting
~~~~~~~~

Use the following command to run the counting step: ::

    $ java -jar Diachromatic.jar count \
        -v prefix.valid_pairs.aligned.bam \
        -d hg19_HinDIII_DigestedGenome.txt


Summarize
~~~~~~~~~

To run the summarize command with the truncate data, run the following command. ::

    $ java -jar target/Diachromatic.jar summarize \
        -o HinD3 \
        -x foo \
        -t HinD3/foo.truncation.stats


This will generate an HTML file called ``HinD3/foo.summary.stats.html``.

