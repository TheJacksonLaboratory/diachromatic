
Tutorial
========

This tutorial shows how to use Diachromatic for processing and quality control of Capture Hi-C reads. Before proceeding with the tutorial, please follow the program setup instructions to build Diachromatic and get bowtie2 as well as the hg19 prebuilt index.


Test dataset
~~~~~~~~~~~~

To get the data, visit this `ftp server <ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/>`_ or use: ::

	wget ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/test_1.fastq
	wget ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/test_2.fastq
	wget ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/hg19_HinDIII_DigestedGenome.txt.gz

Then decompress the digest file: ::

	gunzip hg19_HinDIII_DigestedGenome.txt.gz


Truncation
~~~~~~~~~~

The first step of processing raw FASTQ files with Diachromatic is to recognize and truncate reads with filled-in ligation juctions, which indicate reads that include the junction of the chimeric CHC fragment. This is performed with the truncate subcommand: ::

    $ java -jar Diachromatic.jar truncate \
        -q test_1.fastq \
        -r test_2.fastq \
        -e HinDIII \
        -x prefix \
        -o outdir

.. If a single read is chimeric, it is not possible to map it to one locus, and therefore the 3' portion of the chimeric read is removed ("truncated"), leaving behind the 5' portion of the read that should map to a specific locus. If the 5' sequence is too short to be mapped, the entire read pair is discarded.

.. In practice, only about XXXX percent of the readpairs are truncated.


Mapping
~~~~~~~

The second step of the pipeline is to map the truncated read pairs to the target genome. You also need a file that shows the locations of restriction digests across the genome. This file is included in the test dataset. You can use GOPHER to create probes and the digest file. Diachromatic uses bowtie2 to perform the mapping, and then creates a BAM file containing the valid read pairs.

.. If desired, Diachromatic also outputs BAM files with the discarded (arterfactual or unmappable reads).

Use the following command to run the alignment step: ::

    $ java -jar Diachromatic.jar align \
        -b /usr/bin/bowtie2 \
        -i /path/to/bowtie2index/hg19 \
        -q prefix.truncated_R1.fastq.gz \
        -r prefix.truncated_R2.fastq.gz \
        -d hg19_HinDIII_DigestedGenome.txt \
        -x prefix \
        -o outdir



Note that the bowtie2 index can be downloaded from the
`bowtie 2 site <http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml>`_.  The indices are available as zip archives,
e.g., ``GRCh37.zip``. Once you unzip it, the resulting folder will contain multiple files (GRCh37.1.bt2, GRCh37.3.bt2, GRCh37.rev.1.bt2
GRCh37.2.bt2, GRCh37.4.bt2, GRCh37.rev.2.bt2). You need to pass the path to one of these files without the file suffix.
Assuming the directory is located at ``/some/path/GRCh37``, you would therefore pass ``-i /some/path/GRCh37/CRCh37``.




Counting
~~~~~~~~

Use the following command to run the counting step: ::

    $ java -jar Diachromatic.jar count \
        -v prefix.valid_pairs.aligned.bam \
        -d hg19_HinDIII_DigestedGenome.txt \
        -x prefix \
        -o outdir


Summarize
~~~~~~~~~

To run the summarize command with the truncate data, use the following command. ::

    $ java -jar Diachromatic.jar summarize \
        -t outdir/prefix.truncation.stats.txt \
        -a outdir/prefix.align.stats.txt \
        -c outdir/prefix.count.stats.txt \
        -x prefix \
        -o outdir

This will generate an HTML file called ``outdir/prefix.summary.stats.html``.

The summary results file for the test dataset can also be downloaded from the `ftp server <ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/>`_ or use: ::

	wget ftp://ftp.jax.org/robinp/Diachromatic/test_dataset/test_dataset.summary.stats.html

