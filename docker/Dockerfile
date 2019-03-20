# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make and zlib
RUN apt-get -qqy update \
    && apt-get -qqy install --fix-missing build-essential zlib1g-dev

# Create python environment
RUN conda update -n base conda
RUN conda create --name=olympus python=3.7
RUN conda install -n olympus -y -c bioconda numpy scipy cython natsort boto3 pandas scikit-learn pysam bzip2 pybedtools
ENV CONDA_DEFAULT_ENV olympus
RUN echo "source activate olympus" >> ~/.bashrc
ENV PATH /opt/conda/envs/olympus/bin:$PATH

# Clone Athena repo
RUN git clone https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/athena /opt/athena \
    && cd /opt/athena \
    && git checkout master \
    && cd -

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc \
    && echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]

