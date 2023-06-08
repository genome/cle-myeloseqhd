#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use Cwd qw(abs_path);
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;
use File::Basename;

die "Provide MyeloseqHD output directory" unless @ARGV == 1;
my $dir = $ARGV[0];
#my $dir = '/storage1/fs1/duncavagee/Active/SEQ/MyeloSeqHD/batchdir/CI-733_rerun';
die "$dir is not a valid directory" unless -d $dir;
$dir = abs_path($dir);

my $git_dir = '/storage1/fs1/duncavagee/Active/SEQ/MyeloSeqHD/process/git/cle-myeloseqhd';
my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'MyeloseqHDAnalysis.wdl');
my $json_template = File::Spec->join($git_dir, 'MyeloseqHDAnalysis.json');

my $group  = '/cle/wdl/haloplex';
my $queue  = 'dspencer';
my $docker = 'registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20';

my $user_group = 'compute-gtac-mgi';

my $main_json = File::Spec->join($dir, 'MyeloseqHD.json');
die "$main_json is not valid" unless -s $main_json;

my $key = "RunInfoString"; 
my $result = `grep $key $main_json`;
my ($run_info_str) = $result =~ /\:\s*"(\S+)",/;
#my $run_info_str = '220325_A00118_0477_BH7WJLDSX3,A00118,B,S4,NovaSeqXp,151,10,10,151'; 

my $sample_sheet = File::Spec->join($dir, 'sample_index');
my $fh = IO::File->new($sample_sheet) or die "fail to open $sample_sheet";

my %info;

while (my $line = $fh->getline) {
    chomp $line;
    my @columns = split /\t/, $line;
    $info{$columns[1]} = {
        mrn       => $columns[7],
        all_mrn   => $columns[8],
        accession => $columns[9],
        DOB       => $columns[10],
        sex       => $columns[11],
        exception => $columns[12],
    };
}
$fh->close;

my $ct = 0;
for my $case_dir (glob("$dir/TW*"), glob("$dir/H_*")) {
    my ($case_name) = basename($case_dir) =~ /^(\S+)_[ATCG]{8}$/;
    my $dragen_dir = File::Spec->join($case_dir, 'dragen');
    die "$dragen_dir not existing" unless -d $dragen_dir;
    
    my $bam = File::Spec->join($dragen_dir, $case_name.'_tumor.bam');
    my $bam_index = $bam.'.bai';
    my $vcf = File::Spec->join($dragen_dir, $case_name.'.hard-filtered.vcf.gz');
    my $vcf_index = $vcf.'.tbi';
    
    die "one, some, or all of bam, vcf and their index not valid" unless -s $bam and -s $bam_index and -s $vcf and -s $vcf_index;

    my $inputs = from_json(`cat $json_template`);
    $inputs->{'MyeloseqHDAnalysis.Bam'}            = $bam;
    $inputs->{'MyeloseqHDAnalysis.BamIndex'}       = $bam_index;
    $inputs->{'MyeloseqHDAnalysis.DragenVcf'}      = $vcf;
    $inputs->{'MyeloseqHDAnalysis.DragenVcfIndex'} = $vcf_index;

    $inputs->{'MyeloseqHDAnalysis.MyeloSeqHDRepo'} = $git_dir;

    $inputs->{'MyeloseqHDAnalysis.Name'}           = $case_name;
    $inputs->{'MyeloseqHDAnalysis.OutputDir'}      = $dir;
    $inputs->{'MyeloseqHDAnalysis.SubDir'}         = basename($case_dir);

    $inputs->{'MyeloseqHDAnalysis.mrn'}            = $info{$case_name}->{mrn};
    $inputs->{'MyeloseqHDAnalysis.all_mrn'}        = $info{$case_name}->{all_mrn};
    $inputs->{'MyeloseqHDAnalysis.accession'}      = $info{$case_name}->{accession};
    $inputs->{'MyeloseqHDAnalysis.DOB'}            = $info{$case_name}->{DOB};
    $inputs->{'MyeloseqHDAnalysis.sex'}            = $info{$case_name}->{sex};
    $inputs->{'MyeloseqHDAnalysis.exception'}      = $info{$case_name}->{exception};
    $inputs->{'MyeloseqHDAnalysis.RunInfoString'}  = $run_info_str;

    my $input_json = File::Spec->join($case_dir, 'MyeloseqHDAnalysis.json');
    my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

    $json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
    $json_fh->close;

    my $out_log = File::Spec->join($case_dir, 'out.log');
    my $err_log = File::Spec->join($case_dir, 'err.log');

    my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -R \"select[mem>8000] rusage[mem=8000]\" -M 8000000 -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $input_json $wdl";

    system $cmd;
    #print $cmd."\n";
    print $case_name." submitted\n";
    $ct++;
    #last if $ct == 1;
    sleep 60; #DB upload and query
}
print "All $ct done\n";
