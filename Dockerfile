# Our image is based on Debian bookworm
FROM python:3.11-slim-bookworm

# Set for all apt-get install, must be at the very beginning of the Dockerfile.
ENV DEBIAN_FRONTEND noninteractive

# Get apt dependencies
# - use libatlass instead of liblapack3 libblas3
# - procps is needed for nextflow
RUN apt-get -y update; \
    apt-get -y install \
    build-essential \
    git \
    jq \
    procps \
    && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get -y autoremove

ADD . /app
WORKDIR /app
RUN python setup.py install
# stupid font cache
RUN python -c 'import matplotlib, astropy'

ENTRYPOINT bash

# docker build -t d3vnull0/chips_wrappers:latest . && docker push d3vnull0/chips_wrappers:latest