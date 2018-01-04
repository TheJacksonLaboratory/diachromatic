Diachromatic: Mapping the truncated reads
========================================================================

The next part of the pipeline maps the truncated reads and performs Q/C and filtering.

Mapping procedure
~~~~~~~~~~~~~~~~~
We now align the truncated FASTQ files to the genome. We will start for now with
the test files from the HiCUP distribution, which are human sequences. We will
align to the hg37 genome with bowtie 2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

We download the prebuilt index for the human hg37 genome, H. sapiens, UCSC hg19 (3.5 GB) from the bowtei2 website.
The index needs to be unzipped. ::


    $ unzip hg19.zip
        Archive:  hg19.zip
        inflating: hg19.1.bt2
        inflating: hg19.2.bt2
        inflating: hg19.3.bt2
        inflating: hg19.4.bt2
        inflating: hg19.rev.1.bt2
        inflating: hg19.rev.2.bt2
        inflating: make_hg19.sh

We will call the path to the directory where the index was unpacked /path/to/bowtie2index/




Performing the mapping and Q/C step
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use the following command. ::

    $ java -jar Diachromatic.jar map -b <bowtie2> -i <bowtie2-index> -q <fastq1> -r <fastq2> -d <digest> [-o <outfile>]

The meaning of the options is
    * -b <bowtie2> Path to the bowtie2 executable
    * -i <bowtie2-index> Path to the bowtie2 index for the genome used to map the FASTQ files
    * --q <fastq1> Name and path to the *truncated* "forward" FASTQ file (produced in previous step)
    * --r <fastq2> Name and path to the *truncated* "reverse" FASTQ file (produced in previous step)
    * -d <digest> Path to the digest file produced in the first step
    * [-o <outfile>] This flag is optional and if it is not passed, the default name of ``diachromatic-processed.bam`` will be used.

For instance, xxx. ::

    $ java -jar target/diachromatic-0.0.2.jar map -b /usr/bin/bowtie2 -i btindex/hg19 -q hindIIIhg19chc/test_dataset1.hindIIIhg19.fastq -r hindIIIhg19chc/test_dataset2.hindIIIhg19.fastq -d hg19HindIIIdigest.txtr -o hindIII

Explain output.

The mapped combines mapping and filtering. The steps are
1. Call bowtie2 and map each of the input files separately (as it they were single-end sequences).
2. Then input both of the resulting SAM files. The files are filtered and
mapped accoding to the following criteria
* Discard a read pair if one or more of the reads could not be mapped by bowtei2
* Discard the pair if one of more of the reads was multimapped (bowtie2: XS:i tag exists)
* to do -- filter the read to remove typical HiC artefacts
* todo -- demand that at least one read maps to one of the VPV target regions


