FROM ubuntu:20.04
MAINTAINER David Spencer <dspencer@wustl.edu>

LABEL Image for basic ad-hoc bioinformatic analyses

ARG DEBIAN_FRONTEND=noninteractive

#set timezone to CDT
RUN ln -sf /usr/share/zoneinfo/America/Chicago /etc/localtime
RUN echo "America/Chicago" > /etc/timezone

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    wget \
    g++ \
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libssl-dev \
    libxml2-dev \
    make \
    ncurses-dev \
    nodejs \
    pkg-config \
    unzip \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    bc \
    curl \
    python3.6 \
    python3-pip \
    python3-dev \
    tzdata \
    gawk \
    openssh-client \
    grep \
    evince \
    libnss-sss && \
    apt-get clean all

RUN dpkg-reconfigure --frontend noninteractive tzdata


# Install Python3 packages available through pip 
RUN pip install numpy && \
    pip install cruzdb && \
    pip install cython && \
    pip install pyensembl && \
    pip install pyfaidx && \
    pip install intervaltree_bio && \
    pip install pandas && \
    pip install biopython && \
    pip install scipy && \
    pip install scikit-learn && \
    pip install pysam && \
    pip install statsmodels && \
    pip install pyranges
    
RUN ln -s $(which python3) /usr/local/bin/python

##############
#HTSlib 1.13#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.13.tar.bz2 && \
    cd /tmp/htslib-1.13 && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/ && \
    cp $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix && \
    rm -rf /tmp/htslib-1.13


################
#Samtools 1.13#
################
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar --bzip2 -xf samtools-1.13.tar.bz2 && \
    cd /tmp/samtools-1.13 && \
    ./configure --with-htslib=$HTSLIB_INSTALL_DIR --prefix=$SAMTOOLS_INSTALL_DIR && \
    make && \
    make install && \
    ln -s /opt/samtools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/samtools-1.13


################
#bcftools 1.13#
################
ENV BCFTOOLS_INSTALL_DIR=/opt/bcftools
WORKDIR /tmp
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar --bzip2 -xf bcftools-1.13.tar.bz2 && \
    cd /tmp/bcftools-1.13 && \
    make prefix=$BCFTOOLS_INSTALL_DIR && \
    make prefix=$BCFTOOLS_INSTALL_DIR install && \
    ln -s /opt/bcftools/bin/* /usr/local/bin/ && \
    cd / && \
    rm -rf /tmp/bcftools-1.13

################
#freebayes 1.3.4#
################

WORKDIR /tmp
RUN wget https://github.com/freebayes/freebayes/releases/download/v1.3.4/freebayes-1.3.4-linux-static-AMD64.gz && \
    gunzip freebayes-1.3.4-linux-static-AMD64.gz && \
    cp freebayes-1.3.4-linux-static-AMD64 /usr/local/bin/freebayes && \
    chmod 755 /usr/local/bin/freebayes

################
#tagbam        #
################

RUN apt-get install -y git

WORKDIR /tmp
RUN git clone https://github.com/dhspence/tagbam.git && cd /tmp/tagbam && \
    gcc tagbam.c cgranges.c -o tagbam -lhts -lz && \
    mv tagbam /usr/local/bin && \
    rm -Rf /tmp/tagbam


##############
## bedtools ##
##############

WORKDIR /usr/local
RUN git clone https://github.com/arq5x/bedtools2.git && \
    cd /usr/local/bedtools2 && \
    git checkout v2.27.0 && \
    make && \
    ln -s /usr/local/bedtools2/bin/* /usr/local/bin/

##############
#perl modules#
##############

RUN cpan install Statistics::Basic 
RUN cpan install JSON 

################
#filter script #
################

COPY filterHaloplex.py /usr/local/bin/
COPY CalculateCoverageQC.pl /usr/local/bin/ 
