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
Use the following command: ::

    $ java -jar Diachromatic.jar map -b <bowtie2> -i <bowtie2-index> -q <fastq1> -r <fastq2> -d <digest> [-o <outfile>]

The meaning of the options is:
    * -b <bowtie2> Path to the bowtie2 executable
    * -i <bowtie2-index> Path to the bowtie2 index for the genome used to map the FASTQ files
    * --q <fastq1> Name and path to the *truncated* "forward" FASTQ file (produced in previous step)
    * --r <fastq2> Name and path to the *truncated* "reverse" FASTQ file (produced in previous step)
    * -d <digest> Path to the digest file produced in the first step
    * [-o <outfile>] This flag is optional and if it is not passed, the default name of ``diachromatic-processed.bam`` will be used.
    * [-x] If this is option is used a set, an additional BAM file for rejected pairs will be created. The general tag for rejected reads is ``YY``. See below for tags of individual artifacts.

For instance, the following command will use bowtie2 to map the two FASTQ files of a paired-end run independently (as it they were single-end sequences). Subsequently, the two resulting mappings will be paired, and pairs that show characteristics of known artifacts will be counted and sorted out. Finally, duplicates will be removed. ::

    $ java -jar target/diachromatic-0.0.2.jar map -b /usr/bin/bowtie2 -i btindex/hg19 -q hindIIIhg19chc/test_dataset1.hindIIIhg19.fastq -r hindIIIhg19chc/test_dataset2.hindIIIhg19.fastq -d hg19HindIIIdigest.txtr -o hindIII

Two output files will be produced:
    * ``diachromatic.valid.bam`` contains all uniquely mapped pairs. Known artifacts and duplicated reads are removed. This file can be used for downstream analyses.
    * ``diachromatic.rejected.bam`` contains all pairs that show characteristics of known artifacts:
        * insert too long (Tag: ``TB``)
        * insert too short (Tag: ``TS``)
        * circularized read (Tag: ``SL``)
        * same dangling end (Tag: ``DE``)
        * same internal (Tag: ``SI``)
        * re-ligation (Tag: ``RL``)
        * contiguous (Tag: ``CT``)

Read pairs for which one read cannot be mapped or cannot be mapped uniquely (bowtie2: XS:i tag exists) will be discarded completely. Statistics about the numbers of unmappable reads, multimappable reads, and artifact pairs will be written to the screen.


todo -- demand that at least one read maps to one of the VPV target regions


