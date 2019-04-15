If you dont want or cant install the docker image, there are two other options to install the required dependencies to run the PlasmIdent pipeline.

# Singularity Container

You need a current version (>3.0) of [Singularity](https://www.sylabs.io/singularity/) to download docker images

Singularity containers have advantages over docker and are mostly used in high performance computing contexts. You can use Singularity container by using the singularity profile during launch:

```
nextflow run PlasmIdent/ --input my_input.tsv -profile singularity
```


# Conda environment

The conda environment containing all the dependencies for this pipeline can be installed with the config file `env/PI_env.yml`

You need [miniconda](https://conda.io/miniconda.html) installed on you machine. You can then create the new environment with:

```
conda create -f /env/PI_env.yml 
```

Additionally you need to run the following commands:

1) Fix [bug](https://github.com/bioconda/bioconda-recipes/issues/9830) in Circos:

```
ln -s /opt/conda/envs/PI_env/lib/libwebp.so.6 /opt/conda/envs/PI_env/lib/libwebp.so.7
```

2) Download and install CARD antibiotic resistance database

```
wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate PI_env && rgi load -i card.json
```

You should then be able to run the pipeline without docker by adding the parameter `-profile local`
