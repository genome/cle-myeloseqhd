#!/usr/bin/env python

import sys, os, re, tempfile, csv, pysam, json, binascii, argparse
import pandas as pd
from gtfparse import read_gtf
import pyranges as pr
import numpy as np
from time import gmtime, strftime

#
# Script
#

parser = argparse.ArgumentParser(description='GTF->BED file for exons, with transcript position')
parser.add_argument('-g','--gtf',required=True,help='Ensembl gtf file')
parser.add_argument('-t','--transcripts',required=True,help='List of Ensembl transcripts')

args = parser.parse_args()

gtffile = args.gtf
transcripts = args.transcripts


# read transcripts
trxdf = pd.read_csv(transcripts,sep='\t',names=['gene','transcript'])
 
# read gtf
gtf = read_gtf(gtffile)


for i in trxdf['transcript']:
    df = gtf[(gtf['transcript_id']==i) & (gtf['feature']=='CDS')].copy().reset_index()
    if df['strand'].astype(str).iloc[0] == '-':
        df = df.sort_values('end',ascending=False).reset_index()
        df.at[df.index[-1],'start'] = df.at[df.index[-1],'start']-3
        x = df[['start','end']].apply(lambda row : row['end']-row['start']+1,axis=1).cumsum()
        df['cpos'] = [i+1 for i in [0] + x.tolist()[:-1]]
    else:
        df = df.sort_values('end',ascending=True).reset_index()
        df.at[df.index[-1],'end'] = df.at[df.index[-1],'end']+3
        x = df[['start','end']].apply(lambda row : row['end']-row['start']+1,axis=1).cumsum()
        df['cpos'] = [i+1 for i in [0] + x.tolist()[:-1]]

    df['start'] = df['start']-1
    df = df.sort_values('start',ascending=True)
    
    df[['seqname','start','end','strand','gene_id','gene_name','transcript_id','transcript_name','exon_number','cpos']].to_csv(sys.stdout,sep='\t',header=False, index=False)
