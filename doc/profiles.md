Run profiles
============

Different execution profiles use different parameter settings depending on the computing environment. They are specified during application launch by adding `-profile myProfile` to the command launch.


## Standard

This is the default profile. The pipeline uses docker to download dependencies. It is uses up to 8 cores and 

## Test

```
-profile test
```

You dont need to specify an input file!

This profile runs the whole pipeline on a small sample dataset located in the `data/` folder on only 2 threads. This is mainly for CI testing and debugging purposes but you can use this to check if the installation on your machine was successfull. Uses `docker`

## Local

```
-profile local
```

This profile uses a conda environment with the name `PI_env` installed on your machine. You can create this environment as described [here](alternative_installation.md#conda_environment).

## Singularity

```
-profile singularity
```

This profile uses a singularity container instead of docker. You need to install Singularity


## Binac

```
-profile binac
```

This profile is created to run the pipeline on the BinAc Cluster at Tübingen University using PBS. 

