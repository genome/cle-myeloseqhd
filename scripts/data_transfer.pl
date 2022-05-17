#! /usr/bin/perl

#Copyright (C) 2022 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;
use File::Spec;
use File::Basename;
use File::Copy qw(copy);

umask 002;

unless (@ARGV >= 1) {
    die "MyeloseqHD output directory is required. The string of Sample case dirs that are excluded for transfer are optional"; 
}

my $batch_dir  = $ARGV[0];
my $batch_name = basename $batch_dir;

my $cases_excluded_str;
$cases_excluded_str = $ARGV[1] if @ARGV == 2;
my @excluded_cases = split /,/, $cases_excluded_str;

unless (-d $batch_dir) {
    die "The provided batch dir $batch_dir is not valid";
}

# Run this on scratch1 local
my $staging_dir = "./$batch_name";
if (-d $staging_dir) {
    die "Data transfer staging dir: $staging_dir already existing";
}
else {
    mkdir $staging_dir;
    unless (-d $staging_dir) {
        die "Fail to create $staging_dir";
    }
}

my $QC = File::Spec->join($batch_dir, 'QC_metrics.tsv');
if (-s $QC) {
    copy $QC, $staging_dir;
}
else {
    warn "$QC does not exist or is empty\n";
}

my @files_xfer = qw(
    ampinfo.txt 
    ampcounts.txt 
    cram 
    cram.crai
    haplotectloci.txt 
    haplotect.txt 
    cleaned.vcf.gz
    combined_and_tagged.vcf 
    annotated.vcf.gz
    annotated_filtered.vcf.gz
    report.txt
    report.json
);

opendir(my $dir_h, $batch_dir);
for my $case (readdir $dir_h) {
    next if $case =~ /^(\.|cromwell|dragen|demux|old|test)/;
    next if grep {$_ eq $case} @excluded_cases;
    next unless -d File::Spec->join($batch_dir, $case);

    my ($real_name) = $case =~ /^(\S+lib\d+)_[ATCG]{8}/;
    my $local_tmp_dir = File::Spec->join($staging_dir, $case);
    mkdir $local_tmp_dir;

    for my $file_name (@files_xfer) {
        my $file_path = File::Spec->join($batch_dir, $case, $real_name.'.'.$file_name);
        copy $file_path, $local_tmp_dir;
    }
}
closedir $dir_h;

my $dest_dir = $batch_name.'/';

my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
my $private_key = "/home/$username/.ssh/id_rsa_";
my %info = (
   fdu => {
       suffix => 'pathology',
       user   => 'myloseq_transfer',
   },
   mharriso => { 
       suffix => 'myeloseq',
       user   => 'myloseq_transfer2',
   }
);
$private_key .= $info{$username}->{suffix};

unless (-s $private_key) {
    die "Private key: $private_key is not valid";
}

my $host  = $info{$username}->{user}.'@128.252.17.197';

my $cmd = "/usr/bin/scp -i $private_key -r $staging_dir $host:$dest_dir";
system $cmd;
