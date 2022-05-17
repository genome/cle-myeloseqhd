#!/usr/bin/perl

use POSIX qw/ceil/;

my $covbed = shift;
my $gtf = shift;

`cut -f 9 $covbed | grep ENST | sort -u > /tmp/transcripts.txt`;
`zcat $gtf | awk -v FS=\"\t\" -v OFS=\"\t\" '{ if (!/^#/){ \$1="chr"\$1; } print; }' | grep -w CDS | grep -f /tmp/transcripts.txt > /tmp/transcripts.gtf`;
`awk -v FS=\"\t\" '\$7==\"+\"' /tmp/transcripts.gtf > /tmp/plus_transcripts.gtf`;
`awk -v FS=\"\t\" '\$7==\"-\"' /tmp/transcripts.gtf | sort -k 4,4rn > /tmp/minus_transcripts.gtf`;

my $gene = '';
my $cpos = 0;
open(F,"/tmp/plus_transcripts.gtf");
while(<F>){
    chomp;

    my @l=split("\t",$_);

    /exon_number "(\d+)"/;
    my $ex=$1;
    /gene_name "(\S+)";/;
    my $gn=$1;
    if ($gene ne $gn){
	$cpos=0;
	$gene = $gn;
    }

    foreach my $i ($l[3]..$l[4]){
	$cpos++;
	print join("\t",$l[0],$i-1,$i,$gn,$ex,ceil($cpos / 3)),"\n";
    }

}

open(F,"/tmp/minus_transcripts.gtf");
while(<F>){
    chomp;

    my @l=split("\t",$_);

    /exon_number "(\d+)"/;
    my $ex=$1;
    /gene_name "(\S+)";/;
    my $gn=$1;
    if ($gene ne $gn){
        $cpos=0;
        $gene = $gn;
    }

    foreach my $i (reverse($l[3]..$l[4])){
	$cpos++;
        print join("\t",$l[0],$i-1,$i,$gn,$ex,ceil($cpos / 3)),"\n";
    }

}
