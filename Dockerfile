# Our image is based on Debian bookworm
FROM debian:bookworm

# Set for all apt-get install, must be at the very beginning of the Dockerfile.
ENV DEBIAN_FRONTEND noninteractive

# Get apt dependencies
# - use libatlass instead of liblapack3 libblas3
# - procps is needed for nextflow
RUN apt-get -y update; \
    apt-get -y install \
    build-essential \
    cython3 \
    git \
    ipython3 \
    jq \
    libatlas3-base \
    procps \
    python3-dev \
    python3-ipykernel \
    python3-ipython \
    python3-matplotlib \
    python3-numpy \
    python3-pandas \
    python3-scipy \
    python3-seaborn \
    python3-six \
    python3-astropy \
    wget \
    ; \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    apt-get -y autoremove;

# Bullseye Versions:
# - python3-astropy: 4.2-6
# - python3-h5py: 2.10.0-9
# - python3-numpy: 1:1.19.5-1
# - python3-scipy: 1.6.0-2
# - python3-matplotlib: 3.3.4-1

# Bookworm Versions:
# - python3-astropy: 5.2.1-2
# - python3-numpy: 1.24.2
# - python3-scipy: 1.10.1-2
# - python3-matplotlib: 3.6.3-1


# use python3 as the default python
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 1 \
    && update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1 \
    && update-alternatives --install /usr/bin/ipython ipython /usr/bin/ipython3 1

ADD . /app
WORKDIR /app
RUN python setup.py install

ENTRYPOINT bash