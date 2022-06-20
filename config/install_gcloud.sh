#!/usr/bin/env bash

# Install Google Cloud SDK

apt-get -qqy clean && \
apt-get -qqy update --fix-missing && \
apt-get -qqy dist-upgrade && \
apt-get -qqy install lsb-release && \
pip install -q -U --no-cache-dir \
            crcmod \
            wheel && \
export CLOUD_SDK_REPO="cloud-sdk-$( lsb_release -c -s )" && \
echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get -qqy update && \
apt-get -qqy install --no-install-recommends \
             google-cloud-sdk && \
gcloud config set core/disable_usage_reporting true && \
gcloud config set component_manager/disable_update_check true && \
gcloud config set metrics/environment github_docker_image && \
apt-get -qqy clean && \
rm -rf /tmp/* \
       /var/tmp/* \
       /var/cache/apt/* \
       /var/lib/apt/lists/* \
       /usr/share/man/?? \
       /usr/share/man/??_* && \
gcloud --help
