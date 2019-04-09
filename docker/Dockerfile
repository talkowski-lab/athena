# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make and zlib
RUN apt-get -qqy update \
    && apt-get -qqy install --fix-missing build-essential zlib1g-dev libcurl4-openssl-dev

# Create python environment
# Note: in practice, it's more sane to install packages one at a time
# due to intermittent HTTPS timeouts
RUN conda update -n base conda
RUN conda create --name=olympus python=3.7
RUN conda config --set ssl_verify no
RUN conda install -n olympus -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n olympus -y -c bioconda pysam
RUN conda install -n olympus -y -c bioconda pybedtools
RUN conda install -n olympus -y -c bioconda bzip2
RUN conda install -n olympus -y -c bioconda click
RUN conda install -n olympus -y -c bioconda mysqlclient
RUN conda install -n olympus -y -c bioconda pybigwig
ENV CONDA_DEFAULT_ENV olympus
RUN echo "source activate olympus" >> ~/.bashrc
ENV PATH /opt/conda/envs/olympus/bin:$PATH

# Clone & install Athena
RUN git clone https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/athena /opt/athena \
    && cd /opt/athena \
    && git checkout e52d303cd5156d8c497fb876bcdf471612f80153 \
    && cd -
RUN cd /opt/athena \
    && pip install -e .

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc \
    && echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]

