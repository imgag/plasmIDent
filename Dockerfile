FROM --platform=linux/amd64 condaforge/miniforge3:latest

LABEL maintainer="Caspar Gross <post@caspar.bio>"
LABEL description="Dependencies for the plasmid identification pipeline at github.com/imgag/plasmIDent"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

COPY env/PI_env.yml /tmp/PI_env.yml

RUN apt-get update \
	&& apt-get install --yes --no-install-recommends procps wget bzip2 \
	&& rm -rf /var/lib/apt/lists/*

RUN mamba env create -n PI_env -f /tmp/PI_env.yml \
	&& conda clean --all --force-pkgs-dirs --yes

ENV PATH=/opt/conda/envs/PI_env/bin:$PATH

RUN wget -q -O /tmp/card-data.tar.bz2 https://card.mcmaster.ca/latest/data \
    && mkdir -p /tmp/card_data \
	&& tar -xjf /tmp/card-data.tar.bz2 -C /tmp/card_data \
	&& conda run -n PI_env rgi load --card_json /tmp/card_data/card.json \
	&& rm -rf /tmp/card-data.tar.bz2 /tmp/card_data
