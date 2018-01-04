
Diachromatic: Preprocessing (truncating) raw reads
========================================================================

The next part of the pipeline truncates part of reads that contain the ligation sequence.

Truncation procedure
~~~~~~~~~~~~~~~~~~~~
Valid Hi-C pairs are chimeric reads that require special processing that is distinct from
pipelines for many other genomic analysis techniques such as exome sequencing or ChIP-seq.
The chimeric reads that result from the Hi-C procedure cosist one fragment each from
two different regions of the genome. Capture Hi-C (CHC) is similar except that a capture procedure
is used to enrich regions of interest (usually promoters). Thus, with CHC typical valid
read pairs may comprise a DNA sequence from a promoter and a DNA sequence from an enhancer
that is regulating the promoter.


In most cases, one of the reads of the read pair will map
to a single ligation fragment, and the reverse read will map to another fragment.
However, this is not always true because the Hi-C ligation junction cvan be located within one of the sequenced reads.
This step of the pipeline attempts to address this situation (which could lead to the
read with the Hi-C junction not being mapped during the mapping
step), by deleting sequenced that is downstream of the enzyme recognition site.


For example, if the forward read is entirely contained within one ligation fragment and the reverse read starts in
another fragment but contains the ligation junction and then continues and finishes with part of the fragment of the forward read, then the
truncation step will remove the part of the reverse read that maps
to the first ligation fragment.

Running Diachromatic's Truncation Step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the following command to run the truncation step. ::

    $ java -jar Diachromatic.jar truncate -q test1.fastq -r test2.fastq -e HindIII -s suffix [--outdir <directory>]

The meaning of the options is
    * -q <example1.fastq> Name and path to the "forward" FASTQ file
    * -r <example2.fastq> Name and path to the "reverse" FASTQ file
    * -e <enzyme> The symbol of the restriction enzyme used in the Capture Hi-C experiment
    * -s <suffix> Suffix that will be added to the output files (see below)
    * [--outdir <directory>] Name and path of the output directory. This flag is optional and if it is not passed, the default name of ``results`` will be used.


TODO -- for now I am using the example files from hicup. We should make an example from our
own data!

For instance, the following will truncate FASTQ files from a CHC experiment
that was performed with the enzyme HindIII and place the results into a directory
called hindIIIhg19chc. ::


    $ java -jar target/diachromatic-0.0.2.jar truncate --file1 testdata/test_dataset1.fastq --file2 testdata/test_dataset2.fastq -e HindIII -s hindIIIhg19 --outdir hindIIIhg19chc

If not already present, it creates a directory called ``hindIIIhg19chc``, and writes
two fastq files in it.

    * test_dataset1.hindIIIhg19.fastq
    * test_dataset2.hindIIIhg19.fastq

Note that these files have the same names as the input files except that the suffix ``hindIIIhg19``
was added directory before the final ``fastq`` suffix.
