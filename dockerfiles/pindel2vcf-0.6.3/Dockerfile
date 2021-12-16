FROM ubuntu:xenial
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL \
  version="v0.2.5b8" \
  description="Pindel v0.2.5b8 and samtools 1.9 image for use in Workflows"

RUN apt-get update && apt-get install -y \
    bzip2 \
    libbz2-dev \
    libbz2-1.0 \
    libbz2-dev \
    liblzma-dev \
    liblzma5 \
    liblzma-dev \
    wget \
    make \
    ncurses-dev \
    zlib1g-dev \
    g++

#HTSlib 1.9#
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2 && \
    cd /tmp/htslib-1.9 && \
    ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

#Samtools 1.9#
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar --bzip2 -xf samtools-1.9.tar.bz2 && \
    cd /tmp/samtools-1.9 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.9

# pindel v0.2.5b8
WORKDIR /opt
RUN wget https://github.com/genome/pindel/archive/v0.2.5b8.tar.gz && \
  tar -xzf v0.2.5b8.tar.gz

WORKDIR /opt/pindel-0.2.5b8
RUN ./INSTALL $HTSLIB_INSTALL_DIR

WORKDIR /
RUN ln -s /opt/pindel-0.2.5b8/pindel /usr/local/bin/pindel && \
    ln -s /opt/pindel-0.2.5b8/pindel2vcf /usr/local/bin/pindel2vcf

