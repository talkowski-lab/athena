# Athena: a toolkit for exploring structural variation mutation rates and dosage sensitivity
# Dockerfile
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make and zlib
RUN apt-get -qqy update \
    && apt-get -qqy install --fix-missing \
           apt-utils build-essential zlib1g-dev libcurl4-openssl-dev

# Create python environment
# Note: in practice, it's more sane to install packages one at a time
# due to intermittent HTTPS timeouts
RUN conda update -n base conda
RUN conda create --name=olympus python=3.7
RUN conda config --set ssl_verify no
ENV CONDA_DEFAULT_ENV olympus
RUN echo "source activate olympus" >> ~/.bashrc
ENV PATH /opt/conda/envs/olympus/bin:$PATH
RUN conda install -n olympus -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n olympus -y -c bioconda matplotlib
RUN conda install -n olympus -y -c bioconda samtools=1.9 pysam=0.15.3
RUN conda install -n olympus -y -c bioconda pybedtools bedtools=2.27.0
RUN conda install -n olympus -y -c bioconda bzip2
RUN conda install -n olympus -y -c bioconda click
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing default-libmysqlclient-dev
RUN conda install -n olympus -y -c bioconda mysqlclient
RUN conda install -n olympus -y -c bioconda pybigwig=0.3.17

# Increase max MySQL packet size for large UCSC queries
RUN echo "[mysqld]\nmax_allowed_packet = 2G" > /etc/my.cnf
RUN echo "[client]\nmax_allowed_packet = 2G" >> /etc/my.cnf

# Clone & install Athena
RUN git clone https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/athena /opt/athena \
    && cd /opt/athena \
    && git checkout 2651626a44ba72dd3ab9e630b94d9b06de38a05e \
    && cd -
RUN cd /opt/athena \
    && pip install -e .

# Symlink Athena core reference data to home directory (for GitHub tutorial)
RUN ln -s /opt/athena/data /data

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc \
    && echo 'alias less="zless"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]

