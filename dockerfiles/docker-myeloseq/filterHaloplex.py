from __future__ import division
import sys
import argparse
import pysam
import string
import math
import pandas as pd

from Bio import pairwise2
from scipy.stats import binom
from scipy.stats import chisquare
from scipy.stats import chi2_contingency

def revcomp(seq):
    return seq.translate(maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1]

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('vcffile',help='VCF file')
parser.add_argument('bamfile',help='Consensus BAM file')

parser.add_argument('-r',"--reference",default="/storage1/fs1/dspencer/Active/refdata/GRCh37/all_sequences.fa",
                    help='reference genome fasta file')

parser.add_argument('samplename',help='Sample name')

parser.add_argument('-w',"--window",type=int,default=150,
                    help='window for creating ref and alt sequences')

parser.add_argument("--minreadsperfamily", type=int, default=3,
    help='Minimum number of reads in a family to include it')

parser.add_argument('-n',"--maxreadmismatches",type=int,default=4,
                                        help='Maximum mismatches for a read to be counted')

parser.add_argument('-Q',"--minbasequal",type=int,default=15,
                                        help='Minium base quality to assess variant')

parser.add_argument('-q',"--minmapqual",type=int,default=1,
                                        help='Minium read mapping quality to contribute to variant calling')

parser.add_argument('-m',"--minreads",type=int,default=3,
                    help='minimum reads for defining an amplicon')

parser.add_argument('-v',"--minreads4vaf",type=int,default=3,
                    help='minimum reads in an amplicon for amplicon-based vaf calculation')

parser.add_argument('-a',"--minampnumber",type=int,default=2,
                                        help='Minimum number of covering amplicons with support for variant allele')

parser.add_argument('-p',"--strandpvalue",type=float,default=0.05,
                                        help='Binomial p-value for strand bias')

parser.add_argument('-s',"--minstrandreads",type=int,default=5,
                                        help='Mininum reads on each strand if sufficient strand coverage is present')

parser.add_argument('-d',"--minvaf",type=float,default=0.005,
                                        help='Mininum VAF required')

parser.add_argument('-l',"--lowqualreadbiaspvalue",type=int,default=10,
                                        help='PHRED-scaled low quality read bias p-value cutoff')

# parse arguments
args = parser.parse_args()

####################################
#
# Get general parameters
#
####################################

# sample namne
mysample = args.samplename

# window on either side of a variant for indel annotation 
window = args.window

indelwindow = 50

####################################
#
# Get variant filtering parameters
#
####################################

# minimum calculated VAF
minvaf = args.minvaf

# minimum number of read pairs with a unique insert length to define an amplicon
minreads = args.minreads

# minimum number of reads in an amplicon needed to be eligible to calculate a VAF
minampreadsforvaf = args.minreads4vaf

# minimum number of supporting amplicons
minampnumber = args.minampnumber

# strand bias p value
strandpvalue = args.strandpvalue

# strand read support
minstrandreads = args.minstrandreads

lqrb_pvalue = args.lowqualreadbiaspvalue

####################################
#
# Get read filtering parameters
#
####################################

# min read family size
readfamilysize = args.minreadsperfamily

# maximum number of mismatches in a read to count it (not including read indels)
maxmismatches = args.maxreadmismatches

# min base quality
minqual = args.minbasequal

# minimum mapping quality
minmapqual = args.minmapqual

####################################
#
# Main script
#
####################################

# open vcf file
vcffile = pysam.VariantFile(args.vcffile)
# open bam file
samfile = pysam.AlignmentFile(args.bamfile,"rb")
# open reffasta
fa = pysam.FastaFile(args.reference)

vcffile.header.filters.add("AMPSupport",None,None,'Fails requirement of having >='+str(minampnumber)+' amplicons with support for the variant')
vcffile.header.filters.add("LowQualReadBias",None,None,'PHRED-scaled P-value > 10 of ref and alt alleles in failed vs. passing reads')
vcffile.header.filters.add("LowVAF",None,None,'Calculated VAF <'+str(minvaf))
vcffile.header.filters.add("StrandSupport",None,None,'Variant allele support is lacking on one strand despite adequate coverage (binomial P < '+str(strandpvalue)+' for observing at least '+str(minreads)+' given the coverage on that strand')
#vcffile.header.formats.add("DP", 1, 'Integer', 'Total depth in passing reads')
#vcffile.header.formats.add("AD", 'R', 'Integer', 'Raw/unfiltered number of observation for each allele')
#vcffile.header.formats.add("RO", 1, 'Integer', 'Reference allele observation count in passing reads')
#vcffile.header.formats.add("AO", 1, 'Integer', 'Alternate allele observation count in passing reads')
vcffile.header.formats.add("LQRB", 1, 'String', 'HQ read ref allele count, LQ read ref allele count, HQ read alt allele count, LQ read alt allele count,PHRED-scaled P-value for Low Quality read allele bias')
vcffile.header.formats.add("ST", 1, 'String', 'Counts of plus and minus read counts for alt allele')
vcffile.header.formats.add("TAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position')
vcffile.header.formats.add("SAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position that support the alternate allele')
vcffile.header.formats.add("AMPS", 1, 'String', 'Amplicon string, in the format amplicon id 1,num. reference alleles in amplicon 1, num. alternate alleles in amplicon 1;amplicon id 2,num ref, num alt;...')
vcffile.header.formats.add("CVAF", 1, 'Float', 'Calculated VAF, which is based only on supporting amplicons with with >' + str(minampreadsforvaf) + ' reads')
vcffile.header.add_line("##ampliconparseroptions={'minamplicons':"+str(minampnumber)+"'window':"+str(window)+",'minreads':"+str(minreads)+",'minampreadsforvaf':"+str(minampreadsforvaf)+"}")

hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print("\n".join(hdr) + "\n" + ("\t".join((hd.split("\t"))[0:9])) + "\t" + mysample)

numvars = 0

for vline in vcffile.fetch(reopen=True):

    numvars += 1
    
    # handle multiallelic entries in the VCF (and split them)
    for alt in vline.alts:

        rec = vline
        rec.alts = [ alt ]

        # this is the count of filtered reads supporting ref and alt
        readquals = []
        amplicons = []
        calls = []
        strands = []
        numreads = []
        
        ############################
        #
        # Handle substitutions
        #
        ############################
        if len(rec.ref) == len(rec.alts[0]):

            # for substitutions, use the pileup            
            for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos, max_depth=10000000, multiple_iterators=True):

                # only consider the variant position
                if pileup.pos == rec.pos-1: 
                    
                    for read in pileup.pileups:

                        # count alleles by amplicon
                        if not read.is_del and not read.is_refskip and read.alignment.seq[read.query_position:read.query_position+len(rec.alts[0])] == rec.alts[0]:
                            calls.append('alt')

                        else:
                            calls.append('ref')

    
                        if read.alignment.is_reverse:
                            strands.append('-')
                        else:
                            strands.append('+')

                        # get the amplicon tag
                        if read.alignment.has_tag("XN"):
                            amplicons.append(read.alignment.get_tag("XN"))
                        else:
                            amplicons.append(None)
                            
                        # get read family size
                        numreads.append(int(read.alignment.get_tag("XV")))
                        
                        # skip if more than maxmismatches edit distance for this read or position in indel or 
                        if read.alignment.get_tag("SD") > maxmismatches/read.alignment.query_alignment_length or (not read.is_del and not read.is_refskip and int(read.alignment.query_qualities[read.query_position]) < minqual) or read.alignment.mapping_quality < minmapqual or not read.alignment.has_tag("XN") or not read.alignment.is_proper_pair:
                            readquals.append('fail')

                        else:
                            readquals.append('pass')

        else: # if indel

            # get window around indel
            refseqstart = rec.pos-window-1
            refseqend = rec.pos+window

            # get reference allele sequence
            refseq = fa.fetch(rec.contig,refseqstart,refseqend)

            # indel length
            varlen = max(len(rec.ref) - len(rec.alts[0]),len(rec.alts[0]) - len(rec.ref))

            # get all variants in a window around the indel as a list of tuples => (altseq, is the indel in question)
            altseqs = []
            for var in vcffile.fetch(contig=rec.contig,start=rec.pos-indelwindow,end=rec.pos+indelwindow,reopen=True):
                altseq = '';
                # if deletion
                if len(var.ref) == 1 and len(var.ref) > len(var.alts[0]):
                    altseq = fa.fetch(rec.contig,refseqstart,var.pos-1) + fa.fetch(rec.contig,var.pos+len(var.ref)-1,refseqend+len(var.ref))

                # if complex or insertion or substitution
                else:
                    altseq = fa.fetch(rec.contig,refseqstart,var.pos-1) + var.alts[0] + fa.fetch(rec.contig,var.pos-1+len(var.ref),refseqend)

                # if this is the indel in question add it and make the indicator 1
                if var.pos == rec.pos and var.alts[0] == rec.alts[0]:
                    altseqs.append((altseq,1))

                # otherwise add it and make the indicator 0
                else:
                    altseqs.append((altseq,0))
                    
            # now get all reads that map to the vicinity, including those that are softclipped (+/- 5 bp) and determine the best alignment to ref or alts
            for read in samfile.fetch(rec.contig, rec.pos-1-varlen, rec.pos+varlen,multiple_iterators = True):

                if read.is_unmapped is True:
                    continue
                
                # get read start and end if the whole thing is mapped
                readstart = read.reference_start
                readend = read.reference_end - 1
                if read.cigartuples[0][0] == 4:
                    readstart = readstart - read.cigartuples[0][1]
                
                if read.cigartuples[-1][0] == 4:
                    readend = readend + read.cigartuples[-1][1]

                # if the read overlaps the variant position, or would if the whole thing is aligned, then proceed
                if readstart <= rec.pos-1 and readend >= rec.pos-1:

                    # if there are no indels or softclips, then call it for the reference allele 
                    if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == read.query_length and read.get_cigar_stats()[0][10] < varlen:
                        calls.append('ref')

                    # otherwise, do SW to decide whether the read supports the reference allele or the alt allele
                    else:
                        rdseq = read.query_sequence
            
                        alnref = pairwise2.align.localms(rdseq,refseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)

                        maxscore = 0
                        thisindel = 0
                        for (altseq, isthisindel) in altseqs:                            
                                sa = pairwise2.align.localms(rdseq,altseq, 2, -1, -2, -1,score_only=1,one_alignment_only=1)
                                if sa > maxscore:
                                    maxscore = sa
                                    thisindel = isthisindel
                                    
                        if maxscore > alnref and thisindel == 1:
                            calls.append('alt')

                        else:
                            calls.append('ref')
                            

                    if read.is_reverse:
                        strands.append('-')
                    else:
                        strands.append('+')
    
                    # get the amplicon tag                                                                                                                                    
                    if read.has_tag("XN"):
                        amplicons.append(read.get_tag("XN"))
                    else:
                        amplicons.append(None)

                    # get read family size
                    numreads.append(int(read.get_tag("XV")))
                        
                    # skip if more than maxmismatches edit distance for this read or position in indel or 
                    if read.get_tag("SD") > maxmismatches/read.query_alignment_length or read.mapping_quality < minmapqual or not read.is_proper_pair: #not read.has_tag("XN"):
                        readquals.append('fail')

                    else:
                        readquals.append('pass')


        # now go through the evidence for the variant and print

        # make a dataframe with the calls and evidence from each read
        readdat = pd.DataFrame({
                "calls": calls,
                "amplicons": amplicons,
                "readquals": readquals,
                "numreads": numreads,
                "strands": strands
                })

        # make categoricals
        readdat["calls"] = pd.Categorical(readdat["calls"],categories=["ref","alt"],ordered=True)
        readdat["amplicons"] = pd.Categorical(readdat["amplicons"])
        readdat["strands"] = pd.Categorical(readdat["strands"],categories=["+","-"],ordered=True)

        # get only passing reads with enough read families that have > minreadfamilysize reads
        passing = readdat[(readdat["readquals"]=='pass') & (readdat["numreads"]>=readfamilysize)]

        # num alt and ref allele among passing reads
        ao = passing[(passing["calls"]=='alt')].shape[0]
        ro = passing[(passing["calls"]=='ref')].shape[0]

        # total depth
        dp = passing.shape[0]
        
        # dont even print these
        if ao < minreads:
            continue
        
        # total amplicons
        totalamplicons = passing["amplicons"].cat.remove_unused_categories().value_counts().count()
        
        # get number of amplicons with the alterate allele
        supportingamplicons = passing[(passing["calls"]=='alt')]["amplicons"].cat.remove_unused_categories().cat.categories        

        rawvaf = ao / dp
                
        ampliconcounts = []
        
        # get per amplicon support information
        for a in passing["amplicons"].cat.remove_unused_categories().cat.categories:
            ampliconcounts.append(str(a) + "," + str(passing[(passing["amplicons"]==a) & (passing["calls"]=='ref')].shape[0]) + "," + str(passing[(passing["amplicons"]==a) & (passing["calls"]=='alt')].shape[0]))                

        # chi-squared test on alt allele counts in failed vs. passing reads, if there are ref and alt alleles (ie, some alt reads and not homozygous) 
        twobytwo = [[readdat[(readdat["readquals"]=='pass') & (readdat["calls"]=='ref')].shape[0],readdat[(readdat["readquals"]=='fail') & (readdat["calls"]=='ref')].shape[0]],[readdat[(readdat["readquals"]=='pass') & (readdat["calls"]=='alt')].shape[0],readdat[(readdat["readquals"]=='fail') & (readdat["calls"]=='alt')].shape[0]]]
        failedreadbias = 0
        failedreadbiasstr = str(twobytwo[0][0]) + "," + str(twobytwo[0][1]) + ',' + str(twobytwo[1][0]) + "," + str(twobytwo[1][1]) + ","
        if rec.samples[0]['GT'] != (1,1) and readdat[(readdat["calls"]=='ref')].shape[0] > 0 and readdat[(readdat["calls"]=='alt')].shape[0] > 0 and readdat[(readdat["readquals"]=='pass')].shape[0] > 0 and readdat[(readdat["readquals"]=='fail')].shape[0] > 0:
            failedreadbias = min(99,int(-10 * math.log(chi2_contingency(twobytwo)[1] or 1.258925e-10,10)))

        failedreadbiasstr = failedreadbiasstr + str(failedreadbias)

        # per strand VAFs
        obsvaffwd = 0.0
        obsvafrev= 0.0
        plusaltreads = passing[(passing["calls"]=='alt') & (passing["strands"]=='+')].shape[0]
        minusaltreads = passing[(passing["calls"]=='alt') & (passing["strands"]=='-')].shape[0]
        if passing[(passing["strands"]=='+')].shape[0] > 0:
            obsvaffwd =  plusaltreads / passing[(passing["strands"]=='+')].shape[0]

        if passing[(passing["strands"]=='-')].shape[0] > 0:
            obsvafrev = minusaltreads / passing[(passing["strands"]=='-')].shape[0]

        # do filtering
        rec.filter.clear()
        
        if len(supportingamplicons) < minampnumber:
            rec.filter.add("AMPSupport")
        
        # min vaf filter
        if rawvaf < minvaf:
            rec.filter.add("LowVAF")
            
        if ((passing[(passing["calls"]=='alt') & (passing["strands"]=='+')].shape[0] < minstrandreads and binom.cdf(minstrandreads, passing[(passing["strands"]=='+')].shape[0],obsvaffwd, loc=0) < strandpvalue) or ((passing[(passing["calls"]=='alt') & (passing["strands"]=='-')].shape[0] < minstrandreads and binom.cdf(minstrandreads, passing[(passing["strands"]=='-')].shape[0],obsvafrev, loc=0) < strandpvalue))):
            rec.filter.add("StrandSupport")

        if failedreadbias > lqrb_pvalue:
            rec.filter.add("LowQualReadBias")

        if len(rec.filter.values())==0:
            rec.filter.add("PASS")
            
        mygt = ('.','.')
        for s in rec.samples:
            if rec.samples[s]['GT'] == (1,1) and rawvaf > .99:
                mygt = (1,1)
            else:
                mygt = (0,1)
            
        mysample=0
        gl = (0,0,0)
        if 'GL' in rec.format.keys():
            gl = rec.samples[mysample]['GL']

        rec.samples[mysample].clear()

        rec.samples[mysample]['GT'] = mygt
        rec.samples[mysample]['DP'] = dp
        rec.samples[mysample]['AD'] = (readdat[(readdat["calls"]=='alt')].shape[0],readdat[(readdat["calls"]=='alt')].shape[0])
        rec.samples[mysample]['AO'] = ao
        rec.samples[mysample]['RO'] = ro
        rec.samples[mysample]['GL'] = gl
        rec.samples[mysample]['ST'] = str(plusaltreads) + "," + str(minusaltreads)
        rec.samples[mysample]['LQRB'] = failedreadbiasstr
        rec.samples[mysample]['TAMP'] = int(totalamplicons)
        rec.samples[mysample]['SAMP'] = len(supportingamplicons)
        rec.samples[mysample]['AMPS'] = ';'.join(ampliconcounts)
        rec.samples[mysample]['CVAF'] = rawvaf
    
        print("\t".join(str(rec).rstrip().split("\t")[0:10]))
            
# end vcf fetch

fa.close()
vcffile.close()
samfile.close()
