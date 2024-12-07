FROM ghcr.io/mamba-org/micromamba:2.0.3

ARG USER=testuser
ARG UID=1000
ARG GID=1000

USER root
RUN groupadd -g $GID -o $USER
RUN useradd -m -u $UID -g $GID -o -s /bin/bash $USER

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    git \
    libjansson-dev \
    make \
    wget \
    most \
    && rm -rf /var/lib/apt/lists/*
USER $USER

COPY --chown=$UID:$GID pantas-mamba.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes


ARG MAMBA_DOCKERFILE_ACTIVATE=1

WORKDIR /home/mambauser/
RUN git clone --depth=1 https://github.com/AlgoLab/pantas.git
ENV PATH=$PATH:/home/mambauser/pantas:
WORKDIR /home/mambauser/pantas

RUN wget https://github.com/vgteam/vg/releases/download/v1.56.0/vg
RUN chmod +x vg

RUN git submodule update --init --recursive

WORKDIR  build/deps/sdsl-lite
RUN ./install.sh .

WORKDIR  ../gbwt
RUN make -j4

WORKDIR  ../libbdsg-easy
RUN make -j4

WORKDIR  ../..
RUN g++ -O3 -o annotate annotate.cpp -I./deps/gbwt/include -I./deps/sdsl-lite/include -I./deps/libbdsg-easy/include/ -L./deps/sdsl-lite/lib -L./deps/gbwt/lib/ -L./deps/libbdsg-easy/lib/ -lgbwt -lsdsl -fopenmp -lbdsg -lhandlegraph

WORKDIR  ../..

USER root
RUN mkdir /data && chown -v $UID:$GID /data
ENV DATA_DIR=/data

USER $USER

# VOLUME ["/data"]
#ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "pantas"]
