#!/usr/bin/perl

my $bedtools = '/usr/bin/bedtools';

die "bedtools: $bedtools not found" if !-e $bedtools;

my $exonslop = 2;

my $usage = "$0 <gtf file> <gene name> <Haloplex Amplicon Design File> <hotspot bed file> <coverage bed name> <amplicon bed name>";

# note, hotspot bed file should be in the format: chr start end HOTSPOTQC gene_name exon_number codon_numbers . .

my ($gtf,$genelist,$amplicons,$hotspots,$coverageBed,$ampliconBed) = @ARGV;

die $usage if !-s $gtf or !-s $amplicons or !-s $hotspots or !-s $genelist or !$coverageBed or !$ampliconBed;

my %genes = ();

my %ccdsgenes = ();
my %exons = ();

print STDERR "Loading design file...\n";

open(F,$genelist) || die "cant open file $genelist";
while(<F>){
    chomp;
    $genes{$_} = 1;
    $exons{$_} = 1;
}
close F;

print STDERR "Loading GTF file...\n";

open(F,"gunzip -c $gtf |") || die "Cant open annotation file: $gtf";
while(<F>){
    next if /^#/;
    
    chomp;
    my @l = split("\t",$_);
    
    my %h = ();
    while($l[8]=~/(\S+)\s\"*(\S+?)\"*;/g){
	$h{$1}=$2;
    }
    $h{chr} = $l[0];
    $h{start} = $l[3];
    $h{end} = $l[4];
    $h{strand} = $l[6];

    $h{gene_id} =~ s/\.\d+$//g;
    $h{transcript_id} =~ s/\.\d+$//g;
    
    if ($l[2] eq 'transcript' && defined($genes{$h{gene_name}})){
	
	$genes{$h{gene_name}} = \%h if !defined($genes{$h{gene_name}}{chr}) || (($genes{$h{gene_name}}{end} - $genes{$h{gene_name}}{start}) < ($h{end} - $h{start}));
	
    } elsif ($l[2] eq 'CDS' && defined($exons{$h{gene_name}}) && /appris_principal/){
	$exons{$h{gene_name}} = [] if ref($exons{$h{gene_name}}) ne 'ARRAY';
	push @{$exons{$h{gene_name}}}, \%h;
	
    }
}

print STDERR "Generating coverage QC file...\n";

# generate gene-level bed file
open(B,"| sort -u -k 1,1V -k 2,2n -k 4,4V -k 5,5n | intersectBed -a stdin -b $amplicons -wa | sort -u -k 1,1V -k 2,2n > $coverageBed") or die;
foreach my $g (keys %exons){
    foreach my $e (@{$exons{$g}}){
	print B join("\t",$$e{chr},$$e{start}-1-$exonslop,$$e{end}+$exonslop,
		     $g . '_exon_' . $$e{exon_number},$$e{strand},$g,(($$e{end}+$exonslop) - ($$e{start}-1-$exonslop)),
		     $$e{gene_id},$$e{transcript_id}),"\n";
    }
}

open(HS,$hotspots) || die "cant open hotspot file: $hotspots";
while(<HS>){
    print B $_;
}
close HS;
    
close B;


print STDERR "Generating amplicon bed file...\n";

# now make amplicon file:

`awk -v FS="\t" 'NF==6 && /-/' $amplicons | $bedtools groupby -c 4,5,6 -o first > /tmp/minus.bed`;
`awk -v FS="\t" 'NF==6 && /+/' $amplicons | $bedtools groupby -c 4,5,6 -o first > /tmp/plus.bed`;
`cat /tmp/minus.bed /tmp/plus.bed | sort -k 1,1V -k 2,2n -k 6,6 | bedtools groupby -g 1,2,3,6 -c 4,5 -o first,first | sort -k 1,1V -k 2,2n -k 3,3 | awk -v OFS="\t" '{ print \$1,\$2,\$3,\$5,\$6,\$4; }' > $ampliconBed`;

