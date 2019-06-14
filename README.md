# Diachromatic

## Summary

This package will implements truncation, alignment, artifact filtering and counting of Hi-C read pairs.

## Documentation
Documentation is being prepared for ReadTheDocs. We will not be able to publish the
documentation online until this repository is made public (probably in March.
Therefore, for now generate the documentation with

Documentation is being prepared for ReadTheDocs. The documentation can be generated with:

	$ cd docs
	$ make html

View the documentation as local files: ``diachromatic/docs/_build/html/index.html``


## Program setup

Diachromatic requires Java 8 or higher to run. The source code of Diachromatic can be downloaded from [GitHub](https://github.com/TheJacksonLaboratory/diachromatic). To build the application, clone the repository and create the Java app with maven:

    $ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
    $ cd diachromatic
    $ mvn package

Run Diachromatic to test whether the build process was successful:

    $ java -jar target/Diachromatic.jar

You should see the help message of Diachromatic in the shell.
