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
RUN conda install -n olympus -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n olympus -y -c bioconda pysam bzip2 pybedtools
RUN conda install -n olympus -y -c bioconda click
ENV CONDA_DEFAULT_ENV olympus
RUN echo "source activate olympus" >> ~/.bashrc
ENV PATH /opt/conda/envs/olympus/bin:$PATH

# Clone & install Athena
RUN git clone https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/athena /opt/athena \
    && cd /opt/athena \
    && git checkout e9d202b45bdf58a48ab38a24f44f1a142dcafe5b \
    && cd -
RUN cd /opt/athena \
    && pip install -e .

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc \
    && echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]

