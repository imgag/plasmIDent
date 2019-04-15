################################
#                              #
#         DOCKERFILE           #
#   Hybrid Assembly Pipeline   #
#                              #
################################
FROM continuumio/miniconda3
MAINTAINER Caspar Gross <mail@caspar.one>
LABEL description="contains all the dependencies for plasmid Identification pipeline at github.com/caspargross/plasmidIdentification" 

SHELL ["/bin/bash", "-c"]
# Install procps
RUN apt-get update && apt-get install -y procps

# Install conda envs
ADD env/PI_env.yml /tmp/PI_env.yml
RUN conda env create -f /tmp/PI_env.yml -q && conda clean -a

# Fix bug in circos conda env (https://github.com/bioconda/bioconda-recipes/issues/9830)
RUN ln -s /opt/conda/envs/PI_env/lib/libwebp.so.6 /opt/conda/envs/PI_env/lib/libwebp.so.7

# Download CARD-Antibiotic resistance database
RUN wget -q -O card-data.tar.bz2 https://card.mcmaster.ca/latest/data && tar xfvj card-data.tar.bz2 && bash && source activate PI_env && rgi load -i card.json
