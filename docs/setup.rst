
Program setup
=============

Diachromatic requires Java 8 or higher to run. Diachromatic can be obtained from the Diachromatic `GitHub page <https://github.com/TheJacksonLaboratory/diachromatic>`_. We recommend to download the ``Diachromatic.jar`` file of the latest release on the `release page <https://github.com/TheJacksonLaboratory/diachromatic/releases>`_ of the project.

You can then run the program using this command: ::

    $ java -jar Diachromatic.jar

You should see a help message in the shell.


Building Diachromatic from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To build the application on your own, clone the repository and create the Java app with maven. ::

    $ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
    $ cd diachromatic
    $ mvn package

To test whether the build process was successful, enter the following command: ::

    $ java -jar target/Diachromatic.jar


Preparation for bowtie2
~~~~~~~~~~~~~~~~~~~~~~~

The mapping step of the diachromatic pipeline relies on `bowtie2 <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>`_. If needed, install bowtie2 on your system. For instance, on Debian linux systems bowtie2 can be installed with the following command. ::

	$ sudo apt-get install bowtie2

The prebuilt ``bowtie2`` indices for human hg19 (3.5 GB) and other genome builds can be downloaded from the
`bowtie2 website`_. After downloading the correct archived file to your computer, unpack it with: ::

    $ unzip hg19.zip

In the following pages, we will call the path to the directory where the index was unpacked **/path/to/bowtie2index/**. Substitute this with the actual path on your computer.

.. _bowtie2 website: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

