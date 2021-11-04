FROM ubuntu:xenial
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL \
    description="Image for vep90 GATK3.6 and annotation tools"

RUN apt-get update -y && apt-get install -y \
    apt-utils \
    bioperl \
    bzip2 \
    default-jdk \
    default-jre \
    git \
    libfile-copy-recursive-perl \
    libarchive-extract-perl \
    libarchive-zip-perl \
    libapache-dbi-perl \
    libmodule-build-perl \
    libncurses5-dev \
    make \
    ncurses-dev \
    perl-doc \
    python \
    unzip \
    wget \
    zlib1g-dev

RUN apt-get update -y && apt-get install -y python-pip python-dev build-essential nodejs
RUN pip install --upgrade pip
RUN pip install cyvcf2

RUN ln -s /usr/bin/unzip /bin/unzip

##########
#GATK 3.6#
##########
ENV maven_package_name apache-maven-3.3.9
ENV gatk_dir_name gatk-protected
ENV gatk_version 3.6
RUN cd /tmp/ && wget -q http://mirror.nohup.it/apache/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.zip

# LSF: Comment out the oracle.jrockit.jfr.StringConstantPool.
RUN cd /tmp/ \
    && git clone --recursive https://github.com/broadgsa/gatk-protected.git \
    && cd /tmp/gatk-protected && git checkout tags/${gatk_version} \
    && sed -i 's/^import oracle.jrockit.jfr.StringConstantPool;/\/\/import oracle.jrockit.jfr.StringConstantPool;/' ./public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/varianteval/VariantEval.java \
    && mv /tmp/gatk-protected /opt/${gatk_dir_name}-${gatk_version}
RUN cd /opt/ && unzip /tmp/${maven_package_name}-bin.zip \
    && rm -rf /tmp/${maven_package_name}-bin.zip LICENSE NOTICE README.txt \
    && cd /opt/ \
    && cd /opt/${gatk_dir_name}-${gatk_version} && /opt/${maven_package_name}/bin/mvn verify -P\!queue \
    && mv /opt/${gatk_dir_name}-${gatk_version}/protected/gatk-package-distribution/target/gatk-package-distribution-${gatk_version}.jar /opt/GenomeAnalysisTK.jar \
    && rm -rf /opt/${gatk_dir_name}-${gatk_version} /opt/${maven_package_name}

##############
#HTSlib 1.3.2#
##############
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /tmp/htslib-1.3.2
RUN ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

#######
#tabix#
#######
RUN ln -s $HTSLIB_INSTALL_DIR/bin/tabix /usr/bin/tabix

########
#VEP 90#
########
RUN mkdir /opt/vep/
WORKDIR /opt/vep

RUN git clone https://github.com/Ensembl/ensembl-vep.git
WORKDIR /opt/vep/ensembl-vep
RUN git checkout postreleasefix/90

RUN perl INSTALL.pl --NO_UPDATE

WORKDIR /
RUN ln -s /opt/vep/ensembl-vep/vep /usr/bin/variant_effect_predictor.pl

RUN mkdir -p /opt/lib/perl/VEP/Plugins
COPY Downstream.pm /opt/lib/perl/VEP/Plugins/Downstream.pm
COPY Wildtype.pm /opt/lib/perl/VEP/Plugins/Wildtype.pm

COPY add_annotations_to_table_helper.py /usr/bin/add_annotations_to_table_helper.py

