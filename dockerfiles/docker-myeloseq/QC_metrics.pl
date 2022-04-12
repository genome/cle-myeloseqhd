#! /usr/bin/perl

#Copyright (C) 2022 Feiyu Du <fdu@genome.wustl.edu>
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
use List::Util qw(first);

die "Provide myeloseq output dir" unless @ARGV and @ARGV == 1;
my $dir = $ARGV[0];

unless (-d $dir) {
    die "myeloseq output dir: $dir not existing";
}

my %group1 = (
    TOTAL_BASES        => 'MAPPING/ALIGNING SUMMARY: Total bases',
    TOTAL_READS        => 'MAPPING/ALIGNING SUMMARY: Total input reads',
    MISMATCH_RATE_1    => 'MAPPING/ALIGNING SUMMARY: Mismatched bases R1 (%)',
    MISMATCH_RATE_2    => 'MAPPING/ALIGNING SUMMARY: Mismatched bases R2 (%)',
    PCT_MAPPED_READS   => 'MAPPING/ALIGNING SUMMARY: Mapped reads (%)',
    MEDIAN_INSERT_SIZE => 'MAPPING/ALIGNING SUMMARY: Insert length: median',
    PCT_TARGET_ALIGNED_READS     => 'COVERAGE SUMMARY: Aligned reads in target region (%)',
    AVG_ALIGN_TARGET_COVERAGE    => 'COVERAGE SUMMARY: Average alignment coverage over target region',
    UMI_TARGET_MEAN_FAMILY_DEPTH => 'UMI SUMMARY: On target mean family depth',
    UMI_PCT_READ_FAMILY_LESS_3X  => 'UMI SUMMARY: Read families less than 3x (%)',
    UMI_MEDIAN_READ_FAMILY_SIZE  => 'UMI SUMMARY: Median read family size',
    PCT_LOW_COVERAGE_AMPLICON    => 'AMPLICON SUMMARY: Amplicons with low coverage (%)',
);

my %group2 = (
    VARIANT_COUNTS_TIER_1_2_3 => 'VARIANTCOUNTS: Tier1-3',
    VARIANT_COUNTS_LOW_LEVEL  => 'VARIANTCOUNTS: Low Level',
);

my %group3 = (
    HAPLOTECT_SITES  => 'HAPLOTECT SUMMARY: informativeSites',
    HAPLOTECT_CONTAM => 'HAPLOTECT SUMMARY: contaminationEstimate',
);

my %group4 = (
    FAILED_GENES      => 'FAILED GENES',
    FAILED_GENE_COUNT => 'FAILED GENE COUNT',
    FAILED_EXON_COUNT => 'FAILED EXON COUNT',
);

my %group5 = (
    FAILED_HOTSPOTS => 'FAILED HOTSPOTS',
);

my @headers = ('Case', (sort keys %group1), (sort keys %group2), (sort keys %group3), (sort keys %group4), (sort keys %group5));

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
    my $qc_json_file = $lib_dir."/$case_id.report.json";

    my $json_text = do {
        open(my $json_fh, "<:encoding(UTF-8)", $qc_json_file) or die "fail to open $qc_json_file for $case_name";
        local $/;
        <$json_fh>
    };

    my $json = JSON->new;
    my $data = $json->decode($json_text);

    my @values = ($case_id);

    for my $metric1 (sort keys %group1) {
        my $json_key = $group1{$metric1};
        my ($up_json_key) = $json_key =~ /^(\S+\sSUMMARY):\s/;
        my $value = $data->{QC}->{$up_json_key}->{$json_key};
        if ($metric1 eq 'PCT_LOW_COVERAGE_AMPLICON' or $metric1 eq 'UMI_PCT_READ_FAMILY_LESS_3X') {
           $value = sprintf("%.2f", $value); 
        }
        push @values, $value;
    }
 
    for my $metric2 (sort keys %group2) {
        my ($up_json_key, $json_key) = $group2{$metric2} =~ /^(\S+):\s(.+)$/;
        my $value = $data->{QC}->{$up_json_key}->{$json_key};
        $value = 'NONE' unless $value;
        push @values, $value;
    }

    for my $metric3 (sort keys %group3) {
        my ($up_json_key, $name) = $group3{$metric3} =~ /^(.+):\s(\S+)$/; 
        my @columns = @{$data->{QC}->{$up_json_key}->{columns}};
        my $idx = first {$columns[$_] eq $name} 0..$#columns;
        push @values, $data->{QC}->{$up_json_key}->{data}->[0]->[$idx]; #inside the first array ref
    }  

    for my $metric4 (sort keys %group4) {
        my $json_key = $group4{$metric4};
        my $value = $data->{QC}->{$json_key};
        if ($value eq '') {
            if ($json_key =~ /COUNT/) {
                $value = 0;
            }
            else {
                $value = 'NONE';
            }
        }
        push @values, $value;
    }

    for my $metric5 (sort keys %group5) {
        my $json_key = $group5{$metric5};
        my $value;
        if (%{$data->{QC}->{$json_key}}) {
            my @spots = @{$data->{QC}->{$json_key}->{data}};
            $value = join ',', map{$_->[0].'_'.$_->[2]}@spots;
        }
        else {
            $value = 'NONE';
        }
        push @values, $value;
    }

    $out_fh->print(join "\t", @values);
    $out_fh->print("\n");
    print "$case_id done\n";
}

closedir $dir_h;
