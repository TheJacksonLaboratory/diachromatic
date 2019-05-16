# Diachromatic

Differential analysis of chromatin interactions by capture

## Summary

This package will implement parts of the capture Hi-C pipeline as done in
HiCUP, but will do so in a way that is easily integrated with VPV. The final
part of the pipeline will perform differential analysis of 
promoter reads using a likelihood ratio test.


## Documentation

Documentation is being prepared for ReadTheDocs. The documentation can be generated with:

	$ cd docs
	$ make html

And view the documentation as local files: ``diachromatic/docs/_build/html/index.html``


## Program setup

Diachromatic requires Java 8 or higher to run. The source code of Diachromatic can be downloaded from [GitHub](https://github.com/TheJacksonLaboratory/diachromatic). To build the application, clone the repository and create the Java app with maven:

    $ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
    $ cd diachromatic
    $ mvn package

Run Diachromatic to test whether the build process was successful:

    $ java -jar target/Diachromatic.jar

You should see the help message of Diachromatic in the shell.

