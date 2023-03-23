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
from scipy.stats import fisher_exact

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

parser.add_argument('-u',"--maxnbasesinread",type=int,default=5,
                                        help='Maximum N bases in read to be counted')

parser.add_argument('-Q',"--minbasequal",type=int,default=15,
                                        help='Minium base quality to assess variant')

parser.add_argument('-q',"--minmapqual",type=int,default=1,
                                        help='Minium read mapping quality to contribute to variant calling')

parser.add_argument('-m',"--minreads",type=int,default=5,
                    help='minimum reads supporting a variant to be reported')

parser.add_argument('-v',"--minreads4vaf",type=int,default=3,
                    help='minimum reads in an amplicon for amplicon-based vaf calculation')

parser.add_argument('-a',"--minampnumber",type=int,default=1,
                                        help='Minimum number of covering amplicons with support for variant allele')

parser.add_argument('-p',"--strandpvalue",type=float,default=0.05,
                                        help='Binomial p-value for strand bias')

parser.add_argument('-f',"--fisherstrandbias",type=float,default=30,
                                        help='PHRED-scaled p-value for fisher strand bias filter')

parser.add_argument('-s',"--minstrandreads",type=int,default=5,
                                        help='Mininum reads on each strand if sufficient strand coverage is present')

parser.add_argument('-d',"--minvaf",type=float,default=0.02,
                                        help='Mininum VAF required')

parser.add_argument('-l',"--lowqualreadbiaspvalue",type=int,default=30,
                                        help='PHRED-scaled low quality read bias p-value cutoff')

parser.add_argument("--minhqreads",type=int,default=90,
                                        help='Minimum percent of reads that are high quality at variant position')
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

# strand bias p value
sb_pvalue = args.fisherstrandbias

# strand read support
minstrandreads = args.minstrandreads

lqrb_pvalue = args.lowqualreadbiaspvalue

minhqreads = float(args.minhqreads)/100

####################################
#
# Get read filtering parameters
#
####################################

# min read family size
readfamilysize = args.minreadsperfamily

# maximum number of mismatches in a read to count it (not including read indels)
maxmismatches = args.maxreadmismatches

# maximum number of N basecalls in read to be counted
maxnbasesinread = args.maxnbasesinread

# min base quality
minqual = args.minbasequal

# minimum mapping quality
minmapqual = args.minmapqual

# valid fields (only these will be returned, to minimize clutter):

formats = ('GT','AD','DP','AO','RO','ST','LQRB','TAMP','SAMP','AMPS','VAF')
infotags = ('END','SVLEN','SVTYPE','MyeloSeqHDDB','MyeloSeqHDForceGT')
filters = ('PASS')

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

# remove unneeded formats
for k in vcffile.header.formats.keys():
    if k not in formats:
        vcffile.header.formats[k].remove_header()

# remove unneeded info tags
for k in vcffile.header.info.keys():
    if k not in infotags:
        vcffile.header.info[k].remove_header()

# remove unneeded filter tags
for k in vcffile.header.filters.keys():
    if k not in filters:
        vcffile.header.filters[k].remove_header()

# add info tags (if not already there)
if 'set' not in vcffile.header.info.keys():
    vcffile.header.info.add("set",1,'String','Callers that identified this variant')

# add filters (if not already there)
if 'AMPSupport' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("AMPSupport",None,None,'Fails requirement of having >='+str(minampnumber)+' amplicons with support for the variant')

if 'LowReads' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("LowReads",None,None,'Fails requirement of having >='+str(minreads)+' supporting the variant')

if 'LowQualReadBias' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("LowQualReadBias",None,None,'PHRED-scaled P-value > 10 of non-variant vs variant alleles in failed vs. passing reads')

if 'LowQualReads' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("LowQualReads",None,None,'Failed requirement of having a minimum of 90% high quality reads at the variant position')

if 'LowVAF' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("LowVAF",None,None,'Calculated VAF <'+str(minvaf))

if 'StrandSupport' not in vcffile.header.filters.keys():
    vcffile.header.filters.add("StrandSupport",None,None,'Variant allele support is lacking on one strand despite adequate coverage (binomial P < '+str(strandpvalue)+' for observing at least '+str(minreads)+' given the coverage on that strand')

# add format fields (if not already there)
if 'LQRB' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("LQRB", 1, 'String', 'HQ read non-variant allele count, LQ read non-variant allele count, HQ read variant allele count, LQ read variant allele count,PHRED-scaled P-value for Low Quality read variant allele bias')

if 'ST' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("ST", 1, 'String', 'Counts of plus and minus read counts for alt allele')

if 'TAMP' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("TAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position')

if 'SAMP' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("SAMP", 1, 'Integer', 'Estimated number of Haloplex amplicons at this position that support the alternate allele')

if 'AMPS' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("AMPS", 1, 'String', 'Amplicon string, in the format amplicon id 1,num. reference alleles in amplicon 1, num. alternate alleles in amplicon 1;amplicon id 2,num ref, num alt;...')

if 'VAF' not in vcffile.header.formats.keys():
    vcffile.header.formats.add("VAF", 1, 'Float', 'Variant allele fraction')

if 'AD' not in vcffile.header.formats.keys():
    vcffile.header.formats.add('AD', 'R', 'Integer', 'Allelic depths (counting only informative reads out of the total reads) for the ref and alt alleles in the order listed')

if 'DP' not in vcffile.header.formats.keys():
    vcffile.header.formats.add('DP', 1, 'Integer', 'Approximate read depth (reads with MQ=255 or with bad mates are filtered)')

if 'GT' not in vcffile.header.formats.keys():
    vcffile.header.formats.add('GT', 1, 'String', 'Genotype')

if 'RO' not in vcffile.header.formats.keys():
    vcffile.header.formats.add('RO', 1, 'Integer', 'Reference allele observation count')

if 'AO' not in vcffile.header.formats.keys():
    vcffile.header.formats.add('AO', 1, 'Integer', 'Alternate allele observation count')


# print header
hdr = str(vcffile.header).rstrip().split("\n")
hd = hdr.pop()
print("\n".join(hdr) + "\n" + ("\t".join((hd.split("\t"))[0:9])) + "\t" + mysample)

numvars = 0

for vline in vcffile.fetch(reopen=True):    
    
    numvars += 1

     # first determine the callers that identified the variant
    callers = []
    if 'hotspot' in vline.info.keys():
        callers.append('dragen')

    if 'HOMLEN' in vline.info.keys():
        callers.append('pindel')
   
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
        readids = []

        ############################
        #
        # Handle substitutions
        #
        ############################
        if len(rec.ref) == len(rec.alts[0]):

            # for substitutions, use the pileup
            for pileup in samfile.pileup(rec.contig, rec.pos-1, rec.pos, max_depth=10000000, multiple_iterators=True,ignore_overlaps=False,min_base_quality=0):

                # only consider the variant position
                if pileup.pos == rec.pos-1:

                    for read in pileup.pileups:

                        # get read family size
                        rdsize = int(read.alignment.get_tag("XV"))
                        if rdsize < readfamilysize:
                            continue

                        # check if this read has been counted already.
                        # and only count the read with
                        numreads.append(rdsize)

                        readids.append(read.alignment.query_name)

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
                            amplicons.append("NOTAG")

                        # skip if more than maxmismatches edit distance for this read or position in indel or
                        if read.alignment.get_tag("sd") > maxmismatches/read.alignment.query_alignment_length or read.alignment.seq.count('N') > maxnbasesinread or read.alignment.mapping_quality < minmapqual or not read.alignment.is_proper_pair or (not read.is_del and not read.is_refskip and (int(read.alignment.query_qualities[read.query_position]) < minqual or 'N' in read.alignment.seq[read.query_position:read.query_position+len(rec.alts[0])])):
                            readquals.append('fail')

                        else:
                            readquals.append('pass')

        else: # if indel

            # some pindel calls have N's, which should be skipped all together
            if 'N' in rec.alts[0]:
                continue

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

                altseq = fa.fetch(rec.contig,var.pos-window,var.pos-1) + var.alts[0] + fa.fetch(rec.contig,var.pos-1+len(var.ref),var.pos-1+len(var.ref)+window)

                # if this is the indel in question add it and make the indicator 1
                if var.pos == rec.pos and var.alts[0] == rec.alts[0]:
                    altseqs.append((altseq,1))

                # otherwise add it and make the indicator 0
                else:
                    altseqs.append((altseq,0))

            # now get all reads that map to the vicinity, including those with
            # start and end within 2 * the variant length
            for read in samfile.fetch(rec.contig, rec.pos-1-varlen*2, rec.pos+varlen*2,multiple_iterators = True):
                # get read family size
                rdsize = int(read.get_tag("XV"))
                if read.is_unmapped is True or rdsize < readfamilysize:
                    continue

                # get read start and end if the whole thing is mapped
                readstart = read.reference_start
                readend = read.reference_end - 1

                # if there are softclips, then add varlen slop on either side when determining
                # whether reads overlap the variant position
                if read.cigartuples[0][0] == 4:
                    readstart = readstart - varlen - 1 #read.cigartuples[0][1]

                if read.cigartuples[-1][0] == 4:
                    readend = readend + varlen + 1 #read.cigartuples[-1][1]

                # if the read overlaps the variant position
                if readstart <= rec.pos-1 and readend >= rec.pos-1:

                    readids.append(read.query_name)

                    # if there are no indels or softclips, then call it for the reference allele 
                    if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0 and read.cigartuples[0][1] == read.query_length and read.get_cigar_stats()[0][10] < varlen:
                        calls.append('ref')
                        readquals.append('pass')

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
                            if maxscore >= len(rdseq)*2 *.8:
                                readquals.append('pass')
                            else:
                                readquals.append('fail')

                        else:
                            calls.append('ref')
                            if alnref >= len(rdseq)*2 *.8:
                                readquals.append('pass')
                            else:
                                readquals.append('fail')

                    if read.is_reverse:
                        strands.append('-')
                    else:
                        strands.append('+')

                    # get the amplicon tag or assign notag
                    if read.has_tag("XN"):
                        amplicons.append(read.get_tag("XN"))
                    else:
                        amplicons.append("NOTAG")

                    # get read family size
                    numreads.append(int(read.get_tag("XV")))


        # now go through the evidence for the variant and print

        # make a dataframe with the calls and evidence from each read
        readdat = pd.DataFrame({
                "calls": calls,
                "amplicons": amplicons,
                "readquals": readquals,
                "numreads": numreads,
                "strands": strands,
                "names": readids
                })

        # make categoricals
        readdat["calls"] = pd.Categorical(readdat["calls"],categories=["ref","alt"],ordered=True)
        readdat["amplicons"] = pd.Categorical(readdat["amplicons"])
        readdat["strands"] = pd.Categorical(readdat["strands"],categories=["+","-"],ordered=True)

        # get only passing reads with enough read families that have > minreadfamilysize reads
        passing = readdat[(readdat["readquals"]=='pass') & (readdat["numreads"]>=readfamilysize)]

        # num alt and ref allele among passing reads/read pairs
        ao = passing[(passing["calls"]=='alt')]['names'].drop_duplicates().shape[0]
        ro = passing[(passing["calls"]=='ref')]['names'].drop_duplicates().shape[0]

        # total depth
        dp = passing['names'].drop_duplicates().shape[0]

        # if the variant does exist in the database and there are fewer than minreads or the variant wasnt called then set minreads to 0 because <minreads is below the validated threshold.
        if (ao < minreads or len(callers)==0) and 'MyeloSeqHDDB' in rec.info.keys():
            ao = 0

        # total amplicons
        totalamplicons = passing[['names','amplicons']].drop_duplicates()["amplicons"].cat.remove_unused_categories().value_counts().count()

        # get number of amplicons with the alterate allele
        supportingamplicons = passing[(passing["calls"]=='alt')][['names','amplicons']].drop_duplicates()["amplicons"].cat.remove_unused_categories().cat.categories

        if dp == 0:
            rawvaf = 0
        else:
            rawvaf = ao / dp

        ampliconcounts = []

        # get per amplicon support information
        for a in passing[['names','amplicons']].drop_duplicates()["amplicons"].cat.remove_unused_categories().cat.categories:
            ampliconcounts.append(str(a) + "," + str(passing[(passing["amplicons"]==a) & (passing["calls"]=='ref')][['names','amplicons']].drop_duplicates().shape[0]) + "," + str(passing[(passing["amplicons"]==a) & (passing["calls"]=='alt')][['names','amplicons']].drop_duplicates().shape[0]))

        # chi-squared test on alt allele counts in failed vs. passing reads, if there are ref and alt alleles (ie, some alt reads and not homozygous) and the VAF in failed reads is greater than the VAF in passing reads
        readqual2x2 = [[readdat[(readdat["readquals"]=='pass') & (readdat["calls"]=='ref')].shape[0],readdat[(readdat["readquals"]=='fail') & (readdat["calls"]=='ref')].shape[0]],[readdat[(readdat["readquals"]=='pass') & (readdat["calls"]=='alt')].shape[0],readdat[(readdat["readquals"]=='fail') & (readdat["calls"]=='alt')].shape[0]]]
        passingreadvaf = readqual2x2[1][0]/(readqual2x2[1][0]+readqual2x2[0][0]) if (readqual2x2[1][0]+readqual2x2[0][0]) > 0 else 0
        failingreadvaf = readqual2x2[1][1]/(readqual2x2[1][1]+readqual2x2[0][1]) if (readqual2x2[1][1]+readqual2x2[0][1]) > 0 else 0
        failedreadbias = 0
        failedreadbiasstr = str(readqual2x2[0][0]) + "," + str(readqual2x2[0][1]) + ',' + str(readqual2x2[1][0]) + "," + str(readqual2x2[1][1]) + ","
        if 'GT' in rec.format.keys() and rec.samples[0]['GT'] != (1,1) and readdat[(readdat["calls"]=='ref')].shape[0] > 0 and readdat[(readdat["calls"]=='alt')].shape[0] > 0 and readdat[(readdat["readquals"]=='pass')].shape[0] > 0 and readdat[(readdat["readquals"]=='fail')].shape[0] > 0 and failingreadvaf > passingreadvaf and readqual2x2[1][1]/(readqual2x2[1][1]+readqual2x2[1][0]) > 0.2:
            failedreadbias = min(99,int(-10 * math.log(fisher_exact(readqual2x2)[1] or 1.258925e-10,10)))

        failedreadbiasstr = failedreadbiasstr + str(failedreadbias)

        # per strand VAFs
        obsvaffwd = 0.0
        obsvafrev= 0.0

        plusaltreads = passing[(passing["calls"]=='alt') & (passing["strands"]=='+')].shape[0]
        minusaltreads = passing[(passing["calls"]=='alt') & (passing["strands"]=='-')].shape[0]

        plusrefreads = passing[(passing["calls"]=='ref') & (passing["strands"]=='+')].shape[0]
        minusrefreads = passing[(passing["calls"]=='ref') & (passing["strands"]=='-')].shape[0]

        if passing[(passing["strands"]=='+')].shape[0] > 0:
            obsvaffwd =  plusaltreads / passing[(passing["strands"]=='+')].shape[0]

        if passing[(passing["strands"]=='-')].shape[0] > 0:
            obsvafrev = minusaltreads / passing[(passing["strands"]=='-')].shape[0]

        # Fisher exact SB p-value (phred-scaled)
        # Calculate this only if the ref and alts allele are skewed opposite directions and variant is not homozygous for the alt allele
        sb=0
        sb2x2 = [[plusrefreads,minusrefreads],[plusaltreads,minusaltreads]]
        refstrandskewing = plusrefreads/(plusrefreads+minusrefreads)-.5 if (plusrefreads+minusrefreads) > 0 else 0
        altstrandskewing = plusaltreads/(plusaltreads+minusaltreads)-.5 if (plusaltreads+minusaltreads) > 0 else 0
        if 'GT' in rec.format.keys() and rec.samples[0]['GT'] != (1,1) and rawvaf < .9 and altstrandskewing * refstrandskewing < 0 and readdat[(readdat["calls"]=='ref')].shape[0] > 0 and readdat[(readdat["calls"]=='alt')].shape[0] > 0 and readdat[(readdat["strands"]=='+')].shape[0] > 0 and readdat[(readdat["strands"]=='-')].shape[0] > 0:
            sb = min(min(99,int(-10 * math.log(binom.pmf(plusaltreads, plusaltreads+minusaltreads,plusrefreads/(plusrefreads+minusrefreads)) or 1.258925e-10,10))),min(99,int(-10 * math.log(binom.pmf(plusaltreads, plusaltreads+minusaltreads,minusrefreads/(plusrefreads+minusrefreads)) or 1.258925e-10,10))))

        nrec = vcffile.new_record()

        # if dragen is the only caller then use those counts and the filter flag        
        if len(callers)==1 and callers[0]=='dragen':

            ro = rec.samples[0]['AD'][0]
            ao = rec.samples[0]['AD'][1]
            rawvaf = rec.samples[0]['AF'][0]
            for v in rec.filter.keys():
                nrec.filter.add(v) 

        # if other callers made this call then do custom filtering
        else:
             
            if len(supportingamplicons) < minampnumber and ao > 0:
                nrec.filter.add("AMPSupport")

            # if there are < minstrandreads alt-supporting reads, then calculate the binomial p-value
            # for that observation given the overall VAF and the strand-specific read depth
            if ao > 0 and ((plusaltreads+plusrefreads > 0 and plusaltreads < minstrandreads and binom.cdf(minstrandreads, plusaltreads+plusrefreads,rawvaf, loc=0) < strandpvalue) or (minusaltreads+minusrefreads > 0 and minusaltreads < minstrandreads and binom.cdf(minstrandreads, minusaltreads+minusrefreads,rawvaf, loc=0) < strandpvalue)):
                nrec.filter.add("StrandSupport")

            if ao > 0 and failedreadbias > lqrb_pvalue:
                nrec.filter.add("LowQualReadBias")

            if ao > 0 and readdat.shape[0] > 0:
                if readdat[(readdat["readquals"]=='pass')].shape[0] / readdat.shape[0] < minhqreads:
                    nrec.filter.add("LowQualReads")

        # min vaf filter
        if rawvaf < minvaf:
            nrec.filter.add("LowVAF")

        if len(nrec.filter.values())==0:
            nrec.filter.add("PASS")

        # dont even print variants with <minreads and the variant has no previous records from the database (indicated by the 'MyeloSeqHDDB' info field)
        if ao < minreads and 'MyeloSeqHDDB' not in rec.info.keys() and 'MyeloSeqHDForceGT' not in rec.info.keys():
            continue       

        # otherwise, continue and print variant

        mygt = ('.','.')
        if rawvaf > .98:
            mygt = (1,1)
        else:
            mygt = (0,1)

        mysample=0

        # for sites that require genotyping only
        if 'MyeloSeqHDForceGT' in rec.info.keys() and dp > 100:
            nrec.filter.clear()
            nrec.filter.add("PASS")
            if rawvaf > 0.2 and rawvaf < .98:
                mygt = (0,1)
            elif rawvaf <= .2:
                mygt = (0,0)
            else:
                mygt = (1,1)

        # complete new record
        nrec.chrom = rec.chrom
        nrec.pos = rec.pos
        nrec.alleles = rec.alleles
        nrec.id = rec.id

        # add only specified info tags
        for k in rec.info.keys():
            if k in infotags:
                nrec.info[k] = rec.info.get(k)

        nrec.info['set'] = "-".join(callers) or 'NONE'
        nrec.samples[mysample]['GT'] = mygt
        nrec.samples[mysample]['DP'] = ro+ao
        nrec.samples[mysample]['AD'] = (ro,ao)
        nrec.samples[mysample]['AO'] = ao
        nrec.samples[mysample]['RO'] = ro
        nrec.samples[mysample]['ST'] = str(plusrefreads) + "," + str(minusrefreads) + "," + str(plusaltreads) + "," + str(minusaltreads) + "," + str(sb)
        nrec.samples[mysample]['LQRB'] = failedreadbiasstr
        nrec.samples[mysample]['TAMP'] = int(totalamplicons)
        nrec.samples[mysample]['SAMP'] = len(supportingamplicons)
        nrec.samples[mysample]['AMPS'] = ';'.join(ampliconcounts)
        nrec.samples[mysample]['VAF'] = rawvaf

        # this is a workaround because I cant remove format fields or copy the old record for some reason
        line = str(nrec).rstrip().split("\t")[0:10]
        print("\t".join(line))

# end vcf fetch

fa.close()
vcffile.close()
samfile.close()

