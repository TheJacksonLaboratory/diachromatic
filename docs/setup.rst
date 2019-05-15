
Program setup
=============

Diachromatic requires Java 8 or higher to run. The source code of Diachromatic can be downloaded
from the Diachromatic `GitHub page <https://github.com/TheJacksonLaboratory/diachromatic>`_.

To build the application, clone the repository and create the Java app with maven. ::

    $ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
    $ cd diachromatic
    $ mvn package

To test whether the build process was successful, enter the following command: ::

    $ java -jar target/Diachromatic.jar

You should see a help message in the shell.


Preparation for bowtie2
~~~~~~~~~~~~~~~~~~~~~~~

The mapping step of the diachromatic pipeline relies on `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_. If needed, install bowtie2 on your system. For instance, on Debian linux systems bowtie2 can be installed with the following command. ::

	$ sudo apt-get install bowtie2

The prebuilt ``bowtie2`` indices for human hg19 (3.5 GB) and other genome builds can be downloaded from the
`bowtie2 website`_. After downloading the correct archived file to your computer, unpack it with: ::

    $ unzip hg19.zip

In the following pages, we will call the path to the directory where the index was unpacked **/path/to/bowtie2index/**. Substitute this with the actual path on your computer.

.. _bowtie2 website: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

