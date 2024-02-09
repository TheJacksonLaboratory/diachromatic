# Diachromatic

This package implements truncation, alignment, artifact filtering and counting of Hi-C read pairs.


## Download

Diachromatic requires Java 11 or higher to run. We recommend downloading the latest version from 
the [release page](https://github.com/TheJacksonLaboratory/diachromatic/releases) of the project.

You can run the program using this command:

```shell
java -jar Diachromatic.jar
```

You should see a help message in the shell.


## Documentation

The documentation of Diachromatic is available at https://thejacksonlaboratory.github.io/diachromatic.


## Building Diachromatic

To build the application on your own, clone the repository and build the app with the amazing Maven wrapper:

```shell
git clone https://github.com/TheJacksonLaboratory/diachromatic.git
cd diachromatic
./mvnw package
``` 

To test whether the build process was successful, enter the following command:

```shell
java -jar target/Diachromatic.jar
```
