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

use lib "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/perl5/lib/perl5";
use Spreadsheet::Read;
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;

die "Provide rundir, excel sample spreadsheet, and batch name in order" unless @ARGV == 3;

my ($rundir, $sample_sheet, $batch_name) = @ARGV;
die "$rundir is not valid" unless -d $rundir;
die "$sample_sheet is not valid" unless -s $sample_sheet;

my $dir = '/storage1/fs1/duncavagee/Active/SEQ/MyeloSeqHD';
my $git_dir = File::Spec->join($dir, 'process', 'git', 'cle-myeloseqhd');

my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'MyeloseqHD.wdl');
my $zip  = File::Spec->join($git_dir, 'imports.zip');
my $json_template = File::Spec->join($git_dir, 'MyeloseqHD.json');

my $group  = '/cle/wdl/haloplex';
my $queue  = 'dspencer';
my $docker = 'mgibio/genome_perl_environment:compute1-20';

my $user_group = 'compute-duncavagee';

my $out_dir = File::Spec->join($dir, 'batchdir', $batch_name);
unless (-d $out_dir) {
    unless (mkdir $out_dir) {
        die "Failed to make directory $out_dir";
    }
}

#parsing CoPath dump for MRN, ACCESSION, DOB and Gender
my $copath_dump = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/daily_accession/WML_Daily_Accession_Log_CLE.csv';
my $copath_all_dump = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/daily_accession/WML_Daily_All_Accession_Log_CLE.csv';

my %hash = get_info_hash($copath_dump);
my %all_hash = get_info_hash($copath_all_dump);

#parse sample spreadsheet
my $data = Spreadsheet::Read->new($sample_sheet);
my $sheet = $data->sheet(1);

my $ds_str;
my $si_str;
my $seq_id = 2900000000;

my @cases_excluded;

for my $row ($sheet->rows()) {
    next if $row->[0] =~ /Run|Lane/i;
    unless ($row->[0] =~ /\d+/) {
        die "Lane number is expected, Check sample sheet spreadsheet";
    }
    my ($lane, $flowcell, $name, $index_str, $exception) = @$row;

    if ($exception) {
        $exception =~ s/^\s+//;
        $exception =~ s/\s+$//;
        $exception =~ s/,\s+/,/g;
        $exception =~ s/\s+,/,/g;
        $exception =~ s/\s+/_/g;
    }

    if ($exception and $exception =~ /NOTRANSFER/ and $exception =~ /RESEQ/) {
        die "$name has both NOTRANSFER and RESEQ as exception and only one is allowed.";
    }

    $name =~ s/\s+//g;
    #my ($sample, $mrn, $accession) = $lib =~ /^([A-Z]{4}\-(\d+)\-([A-Z]\d+\-\d+)\-[A-Z0-9]+)\-lib/;
    #my $id = $lib =~ /^H_|Research/ ? 'NONE' : $mrn.'_'.$accession;
    #($sample) = $lib =~ /^(\S+)\-lib/ if $lib =~ /^H_|Research|Clinical|Positive\-Control/;

    my ($lib, $mrn, $accession, $sex, $DOB, $all_MRNs);

    if ($exception and $exception =~ /RESEQ|RESEARCH|NOTRANSFER/) {
        if ($exception =~ /RESEQ/) {
            ($accession) = $name =~ /^(W\d+\-\d+)/;
            die "Fail to get accession number from $name" unless $accession;
            die "Fail to get info from CoPath daily All log for RESEQ: $name" unless exists $all_hash{$accession};
            ($mrn, $sex, $DOB, $all_MRNs) = map{$all_hash{$accession}->{$_}}qw(MRN sex DOB all_MRNs);
            $lib = join '-', 'TWAO', $mrn, $name, 'lib1';
        }
        else {
            $lib = $name;
            unless ($lib =~ /^TW/) {
                my $prefix = 'TWAO-';
                if ($exception =~ /RESEARCH/) {
                    $prefix .= 'Research-' unless $lib =~ /RESEARCH/i;
                }
                $lib = $prefix.$lib;
            }
            $lib = $lib.'-lib1' unless $lib =~ /\-lib/;
            #NOTRANSFER  RESEARCH  They will skip query_DB and upload_DB tasks in WF
            ($mrn, $accession, $sex, $DOB, $all_MRNs) = ('NONE') x 5;
        }
    }
    else {
        ($accession) = $name =~ /^(W\d+\-\d+)/;
        die "Fail to get accession number from $name" unless $accession;
        die "Fail to get info from CoPath daily Active log for $name" unless exists $hash{$accession};
        ($mrn, $sex, $DOB, $all_MRNs) = map{$hash{$accession}->{$_}}qw(MRN sex DOB all_MRNs);
        $lib = join '-', 'TWAO', $mrn, $name, 'lib1';
    }

    my ($index) = $index_str =~ /([ATGC]{8})AT\-AAAAAAAAAA/;

    if ($exception) {
        push @cases_excluded, $lib.'_'.$index if $exception =~ /NOTRANSFER/;
        if ($exception =~ /NOTRANSFER|RESEQ|RESEARCH/) {
            $exception =~ s/,?(NOTRANSFER|RESEQ|RESEARCH),?//;
            $exception = 'NONE' if $exception =~ /^\s*$/;
        }
    }
    else {
        $exception = 'NONE';
    }

    my ($sample) = $lib =~ /^(\S+)\-lib/;

    $ds_str .= join ',', $lane, $lib, $lib, '', $index, '';
    $ds_str .= "\n";
    $si_str .= join "\t", $index, $lib, $seq_id, $flowcell, $lane, $lib, $sample, $mrn, $all_MRNs, $accession, $DOB, $sex, $exception;
    $si_str .= "\n";

    $seq_id++;
}

## Get RunInfoString
my $run_xml = File::Spec->join($rundir, 'RunParameters.xml');
unless (-s $run_xml) {
    die "RunParameters.xml $run_xml is not valid";
}
my $xml_fh = IO::File->new($run_xml) or die "Fail to open $run_xml";
my ($runid, $R1cycle, $R2cycle, $index1cycle, $index2cycle, $fcmode, $wftype, $instr, $side);
my $fctype = 0;

while (my $line = $xml_fh->getline) {
    if ($line =~ /<RunId>(\S+)<\/RunId>/) {
        $runid = $1;
    }
    elsif ($line =~ / ReadName="Read1" Cycles="(\d+)" /) {
        $R1cycle = $1;
    }
    elsif ($line =~ / ReadName="Read2" Cycles="(\d+)" /) {
        $R2cycle = $1;
    }
    elsif ($line =~ / ReadName="Index1" Cycles="(\d+)" /) {
        $index1cycle = $1;
    }
    elsif ($line =~ / ReadName="Index2" Cycles="(\d+)" /) {
        $index2cycle = $1;
    }
    elsif ($line =~ /<InstrumentType>(\S+)<\/InstrumentType>/) {
        $wftype = $1;
    }
    elsif ($line =~ /<InstrumentSerialNumber>(\S+)<\/InstrumentSerialNumber>/) {
        $instr = $1;
    }
    elsif ($line =~ /<Side>(\S+)<\/Side>/) {
        $side = $1;
    }
    elsif ($line =~ /<Type>FlowCell<\/Type>/) {
        $fctype = 1;
    }
    if ($fctype) {
        if ($line =~ /<Mode>(\S+)<\/Mode>/) {
            $fcmode = $1;
            $fctype = 0;
        }
    }
}
$xml_fh->close;

## DRAGEN sample sheet
my $num_N;
if ($index1cycle == 19) {
    $num_N = 11;
}
elsif ($index1cycle == 10) {
    $num_N = 2;
}
else {
    die "Invalid index1cycle $index1cycle";
}

my $dragen_ss  = File::Spec->join($out_dir, 'demux_sample_sheet.csv');
my $ss_fh = IO::File->new(">$dragen_ss") or die "Fail to write to $dragen_ss";
$ss_fh->print("[Settings]\n");
$ss_fh->print("AdapterBehavior,trim\n");
$ss_fh->print("AdapterRead1,AAGATCGGAAGAGCACACGTCTGAACTCC+CAGATCGGAAGAGCACACGTCTGAACTCC+GAGATCGGAAGAGCACACGTCTGAACTCC+TAGATCGGAAGAGCACACGTCTGAACTCC\n");
$ss_fh->print("AdapterRead2,AAAGATCGGAAGAGCGTCGTGTAGGGAAA+CAAGATCGGAAGAGCGTCGTGTAGGGAAA+GAAGATCGGAAGAGCGTCGTGTAGGGAAA+TAAGATCGGAAGAGCGTCGTGTAGGGAAA\n");
$ss_fh->print('OverrideCycles,N1Y150;I8N'.$num_N.";U10;N1Y150\n");
$ss_fh->print("[Data]\n");
$ss_fh->print("Lane,Sample_ID,Sample_Name,Sample_Project,index,index2\n");
$ss_fh->print($ds_str);
$ss_fh->close;

## Sample Index
my $si = File::Spec->join($out_dir, 'sample_index');
my $si_fh = IO::File->new(">$si") or die "Fail to write to $si";
$si_fh->print($si_str);
$si_fh->close;

my $run_info_str = join ',', $runid, $instr, $side, $fcmode, $wftype, $R1cycle, $index1cycle, $index2cycle, $R2cycle; 

## Input JSON
my $inputs = from_json(`cat $json_template`);
$inputs->{'MyeloseqHD.OutputDir'}        = $out_dir;
$inputs->{'MyeloseqHD.IlluminaDir'}      = $rundir;
$inputs->{'MyeloseqHD.SampleSheet'}      = $si;
$inputs->{'MyeloseqHD.DemuxSampleSheet'} = $dragen_ss;
$inputs->{'MyeloseqHD.RunInfoString'}    = $run_info_str;

if (@cases_excluded and $inputs->{'MyeloseqHD.DataTransfer'} eq 'true') {
    $inputs->{'MyeloseqHD.CasesExcluded'} = join ',', @cases_excluded;
}
else {
    delete $inputs->{'MyeloseqHD.CasesExcluded'};
}

my $input_json = File::Spec->join($out_dir, 'MyeloseqHD.json');
my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

$json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
$json_fh->close;

my $out_log = File::Spec->join($out_dir, 'out.log');
my $err_log = File::Spec->join($out_dir, 'err.log');

my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -R \"select[mem>16000] rusage[mem=16000]\" -M 16000000 -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl --imports $zip -i $input_json $wdl";

system $cmd;
#print $cmd."\n";

sub get_info_hash {
    my $file = shift;
    my $dump_fh = IO::File->new($file) or die "fail to open $file for reading";
    my %info_hash;

    while (my $l = $dump_fh->getline) {
        next if $l =~ /^(Accession|,,,,)/;
        next unless $l =~ /MyeloSeq\s/;

        chomp $l;
        $l =~ s/\r//g;

        my @columns = split /,/, $l;
        my $id = $columns[0];
        my $all_MRNs = $columns[8];
        $all_MRNs =~ s/\|/,/g;
        my $sex = $columns[4];

        if ($sex eq 'M') {
            $sex = 'male';
        }
        elsif ($sex eq 'F') {
            $sex = 'female';
        }
        else {
            #Gender U
            warn "WARN: unknown gender $sex found for $id";
        }

        $info_hash{$id} = {
            MRN => $columns[2],
            DOB => $columns[5],
            sex => $sex,
            all_MRNs => $all_MRNs,
        };
    }
    $dump_fh->close;

    return %info_hash;
}
