# Customize the official Debian container to allow Pipeline building of MMG
# David Sherman 2017-02-07

FROM debian

USER root

## Standard build tools
RUN apt-get update && \
    apt-get install -y sudo build-essential git cmake doxygen

## Optional module Scotch
RUN apt-get install -y curl zlib1g-dev
RUN curl -L -O http://gforge.inria.fr/frs/download.php/latestfile/298/scotch_6.0.4.tar.gz && \
    tar xzf scotch_6.0.4.tar.gz && \
    ( cd scotch_6.0.4/src && \
      ln -s Make.inc/Makefile.inc.x86-64_pc_linux2 Makefile.inc && \
      make scotch prefix=/usr && make install ) && \
    rm -rf scotch_6.0.4.tar.gz scotch_6.0.4

## Optional module LinearElasticity
RUN git clone https://github.com/ICStoolbox/Commons.git && \
    mkdir -p Commons/build && \
    ( cd Commons/build && HOME=/usr cmake .. && make && make install ) && \
    git clone https://github.com/ICStoolbox/LinearElasticity.git && \
    mkdir -p LinearElasticity/build && \
    ( cd LinearElasticity/build && HOME=/usr cmake .. && make && make install ) && \
    rm -rf Commons LinearElasticity
