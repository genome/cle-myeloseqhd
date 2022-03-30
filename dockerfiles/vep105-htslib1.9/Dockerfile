FROM ensemblorg/ensembl-vep:release_105.0
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL vep105 and htslib1.9 image

RUN mkdir -p /opt/lib/perl/VEP/Plugins
COPY Downstream.pm /opt/lib/perl/VEP/Plugins/Downstream.pm
