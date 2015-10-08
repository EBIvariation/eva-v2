# European Variation Archive (EVA) v2 [![Build Status](https://travis-ci.org/EBIvariation/eva-v2.svg)](https://travis-ci.org/EBIvariation/eva-v2)

This repository contains work in progress for the next version of the European Variation Archive. If you are looking for the production source code, please check https://github.com/EBIvariation/eva and https://github.com/EBIvariation/eva-pipeline

Its components form the core of the new EVA:

* The pipeline for VCF file processing, based on the Spring Batch framework
* The database access layer

The web services and metadata components will be developed in independent repositories.


## EVA Pipeline
#### Using OpenCGA and Spring batch

The current goal is to allow indexing VCF files into mongoDB.

The approach is to have two different jobs: one for genotyped files, and another for aggregated files.

Both jobs will have four (logical) steps: transformation, loading, statistics and annotation.

The reason for using Spring Batch is tracking job statuses and avoiding waste of computation, as result of repeating just 
the needed steps when something fails, in the more automated way possible.

### Using this tool

You may compile the project with `mvn package` and call the produced jar directly, as `$ java -jar eva-pipeline/target/eva-pipeline-0.1.jar`

We did not implement a custom command line, we are using the `org.springframework.boot.autoconfigure.batch.JobLauncherCommandLineRunner`
class to obtain all the parameters (from `.properties` files or command line). All the parameters you can use are 
documented in `src/main/resources/application.properties`. 

#### Examples

To specify the job to run, there is the `spring.batch.job.names` parameter. You can specify a `.properties` in your 
working directory with this line:

    spring.batch.job.names=variantJob

or in the command line, like this:

    java -jar target/gs-batch-processing-0.1.0.jar --spring.batch.job.names=variantJob

as you see, to set a parameter in the command line, you have to use the name as it appears in the `.properties` after a `--`.




