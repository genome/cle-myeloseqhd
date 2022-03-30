#! /usr/bin/perl

#Copyright (C) 2015 Feiyu Du <fdu@genome.wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use JSON;
use IO::File;
use File::Basename;

die "Provide myeloseq output dir" unless @ARGV and @ARGV == 1;
my $dir = $ARGV[0];

unless (-d $dir) {
    die "myeloseq output dir: $dir not existing";
}

my @case_metrics_list = qw(
    TOTAL_READS 
    MAPPED_READS 
    PERCENT_MAPPED_READS
    ONTARGET_READS 
    PERCENT_ONTARGET_READS
    MEAN_READS_PER_UMI 
    PERCENT_UMI_COVERAGE_OVER_2x
    MEAN_READ_PAIRS_PER_AMPLICON 
    PERCENT_AMPLICON_COVERAGE_0x
    PERCENT_AMPLICON_LOW_COVERAGE
    FAILED_HOTSPOTQC 
    PASSED_HOTSPOTQC
    MEAN_UNIQUE_COVERAGE 
    MEDIAN_UNIQUE_COVERAGE
    PERCENT_UNIQUE_COVERAGE_600x 
    PERCENT_UNIQUE_COVERAGE_2000x
    HAPLOTECT_SITES
    HAPLOTECT_CONTAM
    FAILEDGENES
);

my @headers = ('Case', @case_metrics_list);

my $out_file = $dir.'/QC_metrics.tsv';
my $out_fh = IO::File->new(">$out_file") or die "Failed to write to $out_file";
$out_fh->print(join "\t", @headers);
$out_fh->print("\n");

opendir(my $dir_h, $dir);

for my $case_name (readdir $dir_h) {
    next if $case_name =~ /^(\.|cromwell|dragen|demux|old|test)/;
    my $lib_dir = $dir .'/'.$case_name;
    next unless -d $lib_dir;

    my ($case_id) = $case_name =~ /^(\S+lib\d+)/;
    my $qc_json_file = $lib_dir."/$case_id.qc.json";

    my $json_text = do {
        open(my $json_fh, "<:encoding(UTF-8)", $qc_json_file) or die "fail to open $qc_json_file for $case_name";
        local $/;
        <$json_fh>
    };

    my $json = JSON->new;
    my $data = $json->decode($json_text);

    my @values = ($case_id);

    for my $case_metrics_name (@case_metrics_list) {
        if ($case_metrics_name eq 'FAILED_HOTSPOTQC' or $case_metrics_name eq 'PASSED_HOTSPOTQC') {
            my $info = $data->{$case_metrics_name};
            if (@$info == 0) {
                push @values, 'NA';
            }
            else {
                my $str = join ',', @{$data->{$case_metrics_name}};
                push @values, $str;
            }
        }
        elsif ($case_metrics_name eq 'FAILEDGENES') {
            my $info = $data->{$case_metrics_name};
            if ($info) {
                my @strs;
                for my $gene (sort keys %$info) {
                    push @strs, $gene.'('.sprintf("%.2f",$info->{$gene}).')';
                }
                my $str = join ',', @strs;
                push @values, $str;
            }
            else {
                push @values, 'NA';
            }
        }
        else {
            push @values, $data->{$case_metrics_name};
        }
    }

    $out_fh->print(join "\t", @values);
    $out_fh->print("\n");
    print "$case_id done\n";
}

closedir $dir_h;
