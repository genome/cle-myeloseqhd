#!/usr/bin/perl

use strict;
use Getopt::Long;
  
use File::Basename;
use Statistics::Basic qw(:all);
use JSON;

sub sum {
    my $s = 0;
    map { $s+=$_ } @_;
    $s;
}

sub print_qc {
    my ($id,$dat,$qcDat) = @_;

    my $qcflag = $$qcDat{PASSFLAG};
    my @out = ($id,$dat);
    if (defined($qcDat) && defined($qcDat->{$id})){
	my ($lo,$hi) = @{$qcDat->{$id}};

	if ($hi ne '.' && $lo ne '.'){
	    push @out, '['.join("-",$lo,$hi).']';
	    $qcflag = $qcDat->{FAILFLAG} if $dat > $hi || $dat < $lo;
	} elsif ($hi eq '.' && $lo ne '.'){
	    push @out, "[>$lo]";
	    $qcflag = $qcDat->{FAILFLAG} if $dat < $lo;
	} elsif ($lo eq '.' && $hi ne '.'){
	    push @out, "[<$hi]";
	    $qcflag = $qcDat->{FAILFLAG} if $dat > $hi;
	} else {
	    push @out, "[]";
	}
    } else {
	push @out, "[]";
    }
    return (@out,$qcflag);
}

sub print_gene_qc {
    my ($g,$h,$dat,$qcDat) = @_;

    my $qcflag = $$qcDat{PASSFLAG};
    my @out = ($g);
    foreach my $i (@$h){

	push @out, $dat->{$i};

	if (defined($qcDat) && defined($qcDat->{$i})){
	    my ($lo,$hi) = @{$qcDat->{$i}};
	    if ($hi ne '.' && $lo ne '.'){
	      push @out, '['.join("-",$lo,$hi).']';
	      if($dat->{$i} > $hi || $dat->{$i} < $lo){
		$qcflag = $qcDat->{FAILFLAG};
		$out[$#out] .= '***';
	      }
	    } elsif ($hi eq '.' && $lo ne '.'){
		push @out, "[>$lo]";
		if ($dat->{$i} < $lo){
		  $qcflag = $qcDat->{FAILFLAG};
		  $out[$#out] .= '***';
		}
	    } elsif ($lo eq '.' && $hi ne '.'){
		push @out, "[<$hi]";
		if ($dat->{$i} > $hi){
		  $qcflag = $qcDat->{FAILFLAG};
		  $out[$#out] .= '***';
		}
	    } else {
		push @out, "[]";
	    }
	} else {
	    push @out, "[]";
	} 
    }
    return (@out,$qcflag);
}

my $QCVERSION = basename($0,".pl");

my $qcfail = '***FAILED***';
my $qcpass = '';

my $BEDTOOLS="/usr/local/bin/bedtools";
my $SAMTOOLS="/usr/local/bin/samtools";

my $minBaseQual = 13;
my $minMapQual = 1;

my $DIR = '';
my $REFFASTA='';
my $COVERAGEBED='';
my $QCMETRICS = '';
my $INFOFILE = '';
my $NAME = '';

my $help = '';

my $usage = '';

my $usage = <<END;

  HaloplexQC -r|reference <reference fasta> -d|dir <case directory> -t|targets <target bed file> -q|qcmetrics <qc metrics file> -i|infofile <assay info>

    Optional arguments:

      -t|bedtools [/usr/local/bin/samtools]
      -s|samtools [/usr/local/bin/bedtools]

      -h|help prints help

END

die $usage if $#ARGV < 2;

GetOptions("r|reference=s" => \$REFFASTA,
           "t|targets=s" => \$COVERAGEBED,
	   "d|dir=s" => \$DIR,
	   "n|name=s" => \$NAME,
	   "i|infofile=s" => \$INFOFILE,
	   "q|qc=s" => \$QCMETRICS,
	   "l|bedtools=s" => \$BEDTOOLS,
	   "s|samtools=s" => \$SAMTOOLS,
	   "h|help" => \$help);

die $usage if $help;

die "samtools location not valid: $SAMTOOLS\n" if !-s $SAMTOOLS;
die "bedtools location not valid: $BEDTOOLS\n" if !-s $BEDTOOLS;

die "reference location not valid: $REFFASTA\n" if !-s $REFFASTA;
die "coverage qc bed file location not valid: $COVERAGEBED\n" if !-s $COVERAGEBED;

my $CRAM = `readlink -f $DIR/*.cram`;
chomp $CRAM;
die "consensus bam not valid: $CRAM\n" if !-s $CRAM;

my $VARIANTFILE = `readlink -f $DIR/*.variants_annotated.tsv`;
chomp $VARIANTFILE;
die "variant file not valid: $VARIANTFILE\n" if !-e $VARIANTFILE;

my $HAPLOTECT = `readlink -f $DIR/*.haplotect.txt`;
chomp $HAPLOTECT;
die "variant file not valid: $HAPLOTECT\n" if !-s $HAPLOTECT;

my $HAPLOTECTSITES = `readlink -f $DIR/*.haplotectloci.txt`;
chomp $HAPLOTECTSITES;
die "variant file not valid: $HAPLOTECTSITES\n" if !-e $HAPLOTECTSITES;

my $DRAGENMAPPING = `readlink -f $DIR/dragen/*.mapping_metrics.csv`;
chomp $DRAGENMAPPING;
die "cant find dragen mapping file:: $DRAGENMAPPING\n" if !-s $DRAGENMAPPING;

my $DRAGENUMI = `readlink -f $DIR/dragen/*.umi_metrics.csv`;
chomp $DRAGENUMI;
die "cant find dragen UMI metrics file:: $DRAGENUMI\n" if !-s $DRAGENUMI;

my $DRAGENGC = `readlink -f $DIR/dragen/*.gc_metrics.csv`;
chomp $DRAGENGC;
die "cant find dragen gc metrics file:: $DRAGENGC\n" if !-s $DRAGENGC;

my $DRAGENTARGET = `readlink -f $DIR/dragen/*.target_bed_coverage_metrics.csv`;
chomp $DRAGENTARGET;
die "cant find dragen target mapping file:: $DRAGENTARGET\n" if !-s $DRAGENTARGET;

my $AMPLICONINFO = `readlink -f $DIR/*.ampcounts.txt`;
chomp $AMPLICONINFO;
die "amplicon info file location not valid: $AMPLICONINFO\n" if !-s $AMPLICONINFO;

die "qc file location not valid: $QCMETRICS\n" if !-s $QCMETRICS;
die "info file location not valid: $INFOFILE\n" if !-s $INFOFILE;

# get QC metrics

my %QC = (PASSFLAG => $qcpass, FAILFLAG => $qcfail);
open(Q,$QCMETRICS) || die;
while(<Q>){
    chomp;
    next if /^#/;
    s/\s+$//g;
    my @l = split("\t",$_);
    $QC{$l[0]} = [ $l[1], $l[2] ];
}
close Q;

my @covThresholds = @{$QC{COVERAGE_LEVELS}};
$QC{'PERCENT_TARGET_COVERAGE_'.$covThresholds[0].'x'} = $QC{PERCENT_TARGET_COVERAGE1};
$QC{'PERCENT_TARGET_COVERAGE_'.$covThresholds[1].'x'} = $QC{PERCENT_TARGET_COVERAGE2};
$QC{'PERCENT_UNIQUE_COVERAGE_'.$covThresholds[0].'x'} = $QC{PERCENT_UNIQUE_COVERAGE1};
$QC{'PERCENT_UNIQUE_COVERAGE_'.$covThresholds[1].'x'} = $QC{PERCENT_UNIQUE_COVERAGE2};

my $rg = `$SAMTOOLS view -H $CRAM`;

($rg =~ /ID:(\S+)/) || die "No instrument id found";
my $instrumentid = $1;

($rg =~ /SM:(\S+)/) || die "No sample name found";
my $sample = $1;

($rg =~ /LB:(\S+)/) || die "No library id found";
my $library = $1;

my %out = (sample => $sample,
	   library => $library,
	   instrument_id => $instrumentid,
	   version => $QCVERSION, 
	   lowcov => $covThresholds[0],
	   coverageQCLevels => join(",",@covThresholds),
	   coverage_bed => $COVERAGEBED,
	   reference => $REFFASTA);

#  get dragen mapping metrics
open(DM,"$DRAGENMAPPING") || die "Cannot open dragen mapping metrics file: $DRAGENMAPPING\n";
$out{'MAPPING/ALIGNING SUMMARY'} = {};
while(<DM>){
    chomp;
    my @l = split(",",$_);
    if($l[0] eq 'MAPPING/ALIGNING SUMMARY'){
	$out{'MAPPING/ALIGNING SUMMARY'}{$l[2]} = [$l[3],$l[4] || '.'];
    }
}
close DM;

#  get dragen UMI metrics
$out{'UMI STATISTICS'} = {};
open(DU,"$DRAGENUMI") || die "Cannot open dragen UMI metrics file: $DRAGENUMI\n";
while(<DU>){
    chomp;
    my @l = split(",",$_);
    if($l[2] eq 'Histogram of num supporting fragments' or $l[2] eq 'Histogram of unique UMIs per fragment position'){
	$l[3] =~ s/[\{\}]//g;
	my @c = split /\|/,$l[3];
	$out{'UMI STATISTICS'}{$l[2]} = \@c;	
    } else {
      $out{'UMI STATISTICS'}{$l[2]} =[$l[3],$l[4] || '.'];
    }
}
close DU;
    
#  get dragen GC metrics
open(DG,"$DRAGENGC") || die "Cannot open dragen GC metrics file: $DRAGENGC\n";
$out{'GC METRICS'} = {};
while(<DG>){
    chomp;
    my @l = split(",",$_);
    $out{'GC METRICS'}{$l[2]} = [$l[3],$l[4] || '.'];
}
close DG;

#  get dragen target metrics
open(DT,"$DRAGENTARGET") || die "Cannot open dragen target coverage metrics file: $DRAGENTARGET\n";
$out{'TARGET METRICS'} = {};
while(<DT>){
    chomp;
    my @l = split(",",$_);
    $out{'TARGET METRICS'}{$l[2]} = [$l[3],$l[4] || '.'];
}
close DT;

# get amplicon count info
open(AC,"$AMPLICONINFO") || die "Cannot open amplicon count file: $AMPLICONINFO\n";
$out{'AMPLICON COUNTS'} = {};
while(<AC>){
    chomp;
    my @l = split(/\s/,$_);
    $out{'AMPLICON COUNTS'}{$l[0]} = $l[1];
}
close AC;

my %ampcnt = ();
map { $ampcnt{$_ / 2}++ } values %{$out{'AMPLICON COUNTS'}}; # divide by 2 to get read pairs

$out{TOTAL_READS} = $out{'MAPPING/ALIGNING SUMMARY'}{'Total input reads'}->[0];
$out{MAPPED_READS} = $out{'MAPPING/ALIGNING SUMMARY'}{'Mapped reads'}->[0];
$out{PERCENT_MAPPED_READS} = $out{'MAPPING/ALIGNING SUMMARY'}{'Mapped reads'}->[1];
$out{PROPER_READS} = $out{'MAPPING/ALIGNING SUMMARY'}{'Properly paired reads'}->[0];
$out{PERCENT_PROPER_READS} = $out{'MAPPING/ALIGNING SUMMARY'}{'Properly paired reads'}->[1];
$out{ONTARGET_READS} = $out{'TARGET METRICS'}{'Aligned reads in target region'}->[0];
$out{PERCENT_ONTARGET_READS} = $out{'TARGET METRICS'}{'Aligned reads in target region'}->[1];

$out{MEAN_READS_PER_UMI} = $out{'UMI STATISTICS'}{'Mean family depth'}->[0];

$out{PERCENT_UMI_COVERAGE_OVER_2x} = sprintf("%.1f",(($out{'UMI STATISTICS'}{'Total number of families'}->[0] - ($out{'UMI STATISTICS'}{'Histogram of num supporting fragments'}->[0] + $out{'UMI STATISTICS'}{'Histogram of num supporting fragments'}->[1] + $out{'UMI STATISTICS'}{'Histogram of num supporting fragments'}->[2])) / $out{'UMI STATISTICS'}{'Total number of families'}->[0]) * 100);

$out{MEAN_READ_PAIRS_PER_AMPLICON} = sprintf("%.2f",mean(map { $_ / 2 } values %{$out{'AMPLICON COUNTS'}})); # paired end reads
$out{PERCENT_AMPLICON_COVERAGE_0x} = sprintf("%.1f",$ampcnt{0} / scalar(keys %{$out{'AMPLICON COUNTS'}}) * 100);
$out{PERCENT_AMPLICON_LOW_COVERAGE} = sprintf("%.1f",sum(map { $ampcnt{$_} } 0..($QC{AMPLICON_COVERAGE_LEVELS}->[0]-1)) / scalar(keys %{$out{'AMPLICON COUNTS'}}) * 100);

# get consensus coverage depth for each position in the CoverageQC bed file (that is, the coding portions that will be reported)
my %cov = ();
my %hotspot = ();
my @consensusCoverage = ();
my $totalSize =0;
my %coverageLevels = ();
open(C,"$SAMTOOLS view -T $REFFASTA -b -e '[XV]>2' -f 0x2 $CRAM | $SAMTOOLS depth -d 1000000 -Q $minMapQual -q $minBaseQual -b $COVERAGEBED - | /usr/bin/awk -v OFS=\"\t\" '{ print \$1,\$2-1,\$2,\$3; }' | $BEDTOOLS intersect -a stdin -b $COVERAGEBED -wo |") || die "error running QC script $0";
while(<C>){
    chomp;
    my @l = split("\t",$_);
    if (/HOTSPOTQC/){
	push @{$hotspot{$l[8] . '_codon_' . $l[10]}}, $l[3];
	
    } else {

	$totalSize++;
	push @{$cov{$l[9]}{consensusCoverage}}, $l[3];
	$cov{$l[9]}{geneSize}++;
	map { $cov{$l[9]}{passedCoverage}{$_}++ if $l[3] >= $_ } @covThresholds;

	push @{$cov{$l[9]}{exons}{$l[7]}{consensusCoverage}}, $l[3];
	$cov{$l[9]}{exons}{$l[7]}{exonSize} = $l[6] - $l[5];
	map { $cov{$l[9]}{exons}{$l[7]}{passedCoverage}{$_}++ if $l[3] >= $_ } @covThresholds;

    }
}
close C;

# assay wide coverage metrics
$out{MEAN_UNIQUE_COVERAGE} = sprintf("%.1f",mean(map { @{$cov{$_}{consensusCoverage}} } keys %cov));
$out{MEDIAN_UNIQUE_COVERAGE} = sprintf("%.1f",median(map { @{$cov{$_}{consensusCoverage}} } keys %cov));
foreach my $i (@covThresholds){
    $out{"PERCENT_UNIQUE_COVERAGE_" . $i . "x"} = sprintf("%.1f",sum(map { $cov{$_}{passedCoverage}{$i} } keys %cov) / $totalSize * 100);
}

my @failedexons = ();

# gene level coverage
foreach my $g (sort {$a cmp $b} keys %cov){
    my $con_med = sprintf("%.1f",median(@{$cov{$g}{consensusCoverage}}));
    my $con_mean = sprintf("%.1f",mean(@{$cov{$g}{consensusCoverage}}));

    $out{GENECOV}{$g} = { MEAN_TARGET_COVERAGE => $con_mean, MEDIAN_TARGET_COVERAGE => $con_med,SIZE => $cov{$g}{geneSize} };
    map { $out{GENECOV}{$g}{"PERCENT_TARGET_COVERAGE_" . $_ . "x"} = sprintf("%.1f",$cov{$g}{passedCoverage}{$_} / $cov{$g}{geneSize} * 100) } @covThresholds;

    $out{FAILEDGENES}{$g} = $cov{$g}{passedCoverage}{$covThresholds[0]} / $cov{$g}{geneSize} * 100 if ($cov{$g}{passedCoverage}{$covThresholds[0]} / $cov{$g}{geneSize} * 100 < $QC{'PERCENT_TARGET_COVERAGE_'.$covThresholds[0].'x'}->[0]);
    
    $out{GENECOV}{$g}{FAILED_EXONS} = [];
    foreach my $e (sort keys %{$cov{$g}{exons}}){
	my $con_mean = sprintf("%.1f",mean(@{$cov{$g}{exons}{$e}{consensusCoverage}}));
	$out{GENECOV}{$g}{exons}{$e} = { MEAN_TARGET_COVERAGE => $con_mean,
					 MEDIAN_TARGET_COVERAGE => $con_med,
					 SIZE => $cov{$g}{exons}{$e}{exonSize} };
	map { $out{GENECOV}{$g}{exons}{$e}{"PERCENT_TARGET_COVERAGE_" . $_ . "x"} = sprintf("%.1f",$cov{$g}{exons}{$e}{passedCoverage}{$_} / $cov{$g}{exons}{$e}{exonSize} * 100) } @covThresholds;

	push @{$out{GENECOV}{$g}{FAILED_EXONS}}, $e if $cov{$g}{exons}{$e}{passedCoverage}{$covThresholds[0]} / $cov{$g}{exons}{$e}{exonSize} * 100 < $QC{'PERCENT_TARGET_COVERAGE_'.$covThresholds[0].'x'}->[0];
    }
    $out{GENECOV}{$g}{FAILED_EXON_COUNT} = scalar @{$out{GENECOV}{$g}{FAILED_EXONS}};
}

# hotspot QC
my @passed = ();
my @failed = ();
foreach my $i (sort keys %hotspot){
    if ((sort {$a<=>$b} @{$hotspot{$i}})[0] < $covThresholds[0]){
	push @failed, $i;
    } else {
	push @passed, $i;
    }
}

$out{FAILED_HOTSPOTQC} = \@failed;
$out{PASSED_HOTSPOTQC} = \@passed;
$out{FAILED_EXONS} = [ map { @{$out{GENECOV}{$_}{FAILED_EXONS}} } sort keys %{$out{GENECOV}} ];

# get variant information 
my @data;

format qc_format =
@<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     @<<<<<<<<<<<<<< @<<<<<<<<<<<<<<<     @<<<<<<<<<<<<<
@data
.

format gene_qc_format_header =
@<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  @<<<<<<<<<<<<<<<<<<<<<<<<<<<
@data
.

format gene_qc_format =
@<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<<<<< @<<<<<<<<<<< @<<<<<<<<<<<<<<<<  @<<<<<<<<<<< @<<<<<<<<<<<<<< @<<<<<<<<<<
@data
.

$out{timestamp} = scalar localtime();


# print QC report
open(F,">$NAME.qc.txt") || die;
select(F);
$~="qc_format";

print F "SAMPLE ID: $out{sample}\tLIBRARY ID: $out{library}\tINSTRUMENT ID: $out{instrument_id}\tDATE REPORTED: ", $out{timestamp}, "\n\n";

print F "SEQUENCE DATA QC:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } qw(TOTAL_READS MAPPED_READS PERCENT_MAPPED_READS PROPER_READS PERCENT_PROPER_READS ONTARGET_READS PERCENT_ONTARGET_READS);
print F "\n";

print F "LIBRARY QC:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } qw(MEAN_READS_PER_UMI PERCENT_UMI_COVERAGE_OVER_2x);
print F "\n";

map { @data = print_qc($_,$out{$_},\%QC); write; } qw(MEAN_READ_PAIRS_PER_AMPLICON PERCENT_AMPLICON_LOW_COVERAGE);
print F "\n";

print F "HOTSPOT QC:\n\n";
@data = ("FAILED_HOTSPOTQC",scalar @{$out{FAILED_HOTSPOTQC}},"[0]", scalar @{$out{FAILED_HOTSPOTQC}} > $QC{FAILED_HOTSPOTQC}[0] ? $qcfail : $qcpass);
write;

@data = ("PASSED_HOTSPOTQC",scalar @{$out{PASSED_HOTSPOTQC}},"[".(scalar keys %hotspot). "]",scalar @{$out{PASSED_HOTSPOTQC}} < $QC{PASSED_HOTSPOTQC}[0] ? $qcfail : $qcpass);
write;
print F "\n";

if (scalar @{$out{FAILED_HOTSPOTQC}} > 0){
  print F "**FAILED_HOTSPOTS**\t ". join(",",@{$out{FAILED_HOTSPOTQC}}),"\n\n";
}

print F "ASSAY-WIDE COVERAGE METRICS:\n\n";
map { @data = print_qc($_,$out{$_},\%QC); write; } (qw(MEAN_UNIQUE_COVERAGE),map { "PERCENT_UNIQUE_COVERAGE_" . $_ . "x" } @covThresholds);
print F "\n";

print F "GENE/TARGET LEVEL COVERAGE QC:\n\n";

my @headers = (qw(MEAN_TARGET_COVERAGE),(map { "PERCENT_TARGET_COVERAGE_" . $_ . "x" } @covThresholds),"FAILED_EXON_COUNT");
$~="gene_qc_format_header";
@data = ("GENE",@headers);
write;

$~="gene_qc_format";
foreach my $g (sort keys %{$out{GENECOV}}){
    @data = print_gene_qc($g,\@headers,$out{GENECOV}{$g},\%QC);
    write;
}
print F "\n";

print F "FAILED EXONS:\n\n";
@headers = (qw(MEAN_TARGET_COVERAGE),(map { "PERCENT_TARGET_COVERAGE_" . $_ . "x" } @covThresholds));
$~="gene_qc_format_header";
@data = ("EXON",@headers,"");
write;

$~="gene_qc_format";
foreach my $g (sort keys %{$out{GENECOV}}){
    foreach my $e (sort @{$out{GENECOV}{$g}{FAILED_EXONS}}){
	@data = (print_gene_qc($e,\@headers,$out{GENECOV}{$g}{exons}{$e},\%QC),"","","","","");
        write;
    }
}
print F "\n";


my %aa3to1 = qw(Ala A Arg R Asn N Asp D Asx B Cys C Glu E Gln Q Glx Z Gly G His H Ile I Leu L Lys K Met M Phe F Pro P Ser S Thr T Trp W Tyr Y Val V Xxx X Ter *);

my %vars = ("TIER 1-2" => [], "TIER 3" => [], "TIER 4" => [], "FILTERED" => []);
open(VF,$VARIANTFILE) || die;
$_ = <VF>;
chomp;
s/\s+$//;
s/\S+\.CVAF/VAF/;
s/\S+\.TAMP/AMPLICONS/;
s/\S+\.SAMP/SUPPORTING_AMPLICONS/;
s/\S+\.RO/REFERENCE_READS/;
s/\S+\.AO/VARIANT_READS/;

@headers = split("\t",$_);

while(<VF>){
    chomp;
    s/\s+$//;
    my @F = split("\t",$_);
    my %F = map { $headers[$_] => $F[$_] } 0..$#F;

    $F{Consequence} = (split("&",$F{Consequence}))[0];
    $F{Consequence} =~ /^(\S+?)_/;
    $F{EXON} =~ s/\/\d+//g;
    $F{INTRON} =~ s/\/\d+//g;

    my $HGVSpShort = $F{HGVSp};
    $HGVSpShort =~ s/\S+:p\.(\S+?)/\1/g;

    while( $HGVSpShort and my ( $find, $replace ) = each %aa3to1 ) {
        eval "\$HGVSpShort =~ s{$find}{$replace}g";
    }
    $F{HGVSp} = $HGVSpShort;

    my $cat = ($F{MYELOSEQ_MDS_AC}>0 || $F{MYELOSEQ_TCGA_AC}>0) ? "TIER 1-2" : "TIER 3";
    $cat = "TIER 4" if $F{Consequence}=~/synonymous/;

    next if $cat eq 'TIER 4' or $F{Consequence} eq 'synonymous_variant';

    # special case to handle FLT3 ITD
    if ($F{SYMBOL} eq 'FLT3' && $F{EXON} == 14 && length($F{ALT})>length($F{REF}) && length($F{ALT}) > 6){
        $cat = "TIER 1-2" ;
    }

    $cat = "FILTERED" if $F{FILTER} ne 'PASS';
    $cat = "LOWVAF" if $F{FILTER} eq 'LowVAF';

    $F{MAX_AF} = $F{MAX_AF} ? sprintf("%.3f\%",$F{MAX_AF} * 100) : 'none';

    push @{$vars{$cat}}, [$cat,$F{FILTER},$F{SYMBOL},$F{CHROM},$F{POS},$F{REF},$F{ALT},$F{Consequence},$F{HGVSp},$F{HGVSc},$F{EXON} || 'NA',$F{INTRON} || 'NA',$F{MAX_AF},$F{REFERENCE_READS}+$F{VARIANT_READS},$F{VARIANT_READS},sprintf("%.2f\%",$F{VAF} * 100),$F{AMPLICONS},$F{SUPPORTING_AMPLICONS}];

}
close VF;

$out{VARIANTS} = \%vars;

print F "TIER 1-3 Variant detail:\n\n";
if (scalar @{$vars{"TIER 1-2"}} > 0 || scalar @{$vars{"TIER 3"}} > 0){
   print F uc(join("\t", qw(category filter gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))), "\n";
   print F join("\n", map { join("\t", @{$_}) } (@{$vars{"TIER 1-2"}},@{$vars{"TIER 3"}})), "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "LowVAF Variants:\n\n";
if (scalar @{$vars{"LOWVAF"}} > 0){
   print F uc(join("\t", qw(category filter gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))), "\n";
   print F join("\n", map { join("\t", @{$_}) } @{$vars{"LOWVAF"}}), "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "Filtered Variants:\n\n";
if (scalar @{$vars{"FILTERED"}} > 0){
   print F uc(join("\t", qw(category filter gene chrom position ref alt consequence p_syntax c_syntax exon intron pop_freq coverage variant_reads vaf amplicons supporting_amplicons))),"\n";
   print F join("\n", map { join("\t", @{$_}) } @{$vars{"FILTERED"}}) . "\n\n";
} else {
 print F "NO VARIANTS\n\n";
}

print F "CONTAMINATION ESTIMATE\n\n";
open(H,$HAPLOTECT) || die "Cant open haplotect file: $HAPLOTECT";
$_ = <H>;
print F $_;
$_ = <H>;
chomp;
close H;
my @h = split("\t",$_);
print F join("\t",(@h,($h[2] >= $QC{HAPLOTECT_SITES} && $h[6] > $QC{HAPLOTECT_CONTAM} ? $qcfail : $qcpass))),"\n";
$out{HAPLOTECT_SITES} = $h[2];
$out{HAPLOTECT_CONTAM} = $h[6];
$out{HAPLOTECT_CONTAMSNP} = $h[7];

print F "\n\n";

print F "CONTAMINATING HAPLOTYPE INFORMATION\n\n";
open(H,$HAPLOTECTSITES) || die "Cant open haplotect sites file: $HAPLOTECTSITES";
while(<H>){
  print F $_;
}
close H;

print F "\n\n";

open(R,$INFOFILE) || die "Cant open report info file: $INFOFILE";
while(<R>){
    print F $_;
}
close R;
print F "\n";

print F "#qc version: $out{version} mincov value: $out{lowcov}\n";
print F "#target bed file: $out{target_bed}\n";

close F;


# print QC json file
open(F,">$NAME.qc.json") || die;
print F encode_json \%out, "\n";
close F;
