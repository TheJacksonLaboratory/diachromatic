# Diachromatic

This package implements truncation, alignment, artifact filtering and counting of Hi-C read pairs. The documentation of Diachromatic is available on `ReadTheDocs <https://diachromatic.readthedocs.io/en/latest/>`_.


## Download

Diachromatic requires Java 8 or higher to run. Diachromatic can be obtained from the Diachromatic `GitHub page <https://github.com/TheJacksonLaboratory/diachromatic>`_. We recommend to download the ``Diachromatic.jar`` file of the latest release on the `release page <https://github.com/TheJacksonLaboratory/diachromatic/releases>`_ of the project.

You can run the program using this command:

	$ java -jar Diachromatic.jar

You should see a help message in the shell.


## Building Diachromatic

To build the application on your own, clone the repository and create the Java app with maven:

	$ git clone https://github.com/TheJacksonLaboratory/diachromatic.git
	$ cd diachromatic
	$ mvn package

To test whether the build process was successful, enter the following command:

	$ java -jar target/Diachromatic.jar


## Documentation

The documentation of Diachromatic is available on `ReadTheDocs <https://diachromatic.readthedocs.io/en/latest/>`_. You can also generate the documentation yourself with:

	$ cd docs
	$ make html

View the documentation as local files: ``diachromatic/docs/_build/html/index.html``

