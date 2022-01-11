FROM mambaorg/micromamba:latest
USER root
RUN apt-get update && apt-get install -y git acl libxml2-dev build-essential libcurl4-openssl-dev pkg-config m4 libtool automake autoconf libjson-c-dev libgl1-mesa-dev libgl-dev libc++-11-dev libc++abi-11-dev
RUN setfacl -Rm u:micromamba:rwx /opt
USER micromamba
ENV BASE=/home/micromamba
WORKDIR $BASE
RUN git clone https://github.com/mittinatten/freesasa
WORKDIR $BASE/freesasa
RUN git submodule init && git submodule update
RUN autoreconf -i
RUN ./configure
RUN make
USER root
RUN make install
USER micromamba
ADD ./ $BASE/haddock-tools
WORKDIR $BASE/haddock-tools
RUN micromamba env create -y -f $BASE/haddock-tools/environment.yaml
USER root
RUN make
USER micromamba
ENV ENV_NAME="haddock-tools"
ENV PATH="$PATH:$BASE/haddock-tools:$BASE/freesasa"