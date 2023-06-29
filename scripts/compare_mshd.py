import json, os, glob, sys
import pandas as pd
import numpy as np
import argparse

def search_files(directory, pattern):
    file_paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if glob.fnmatch.fnmatch(file_path, pattern):
                file_paths.append(file_path)
    return file_paths

def compare_cases(row):

    out = pd.DataFrame(columns=['case','reference','query','reference_reportable_variants','query_reportable_variants','reference_unique_reportable_variants','query_unique_reportable_variants','TP','FN','FP'])
  
    refj = pd.read_json(row['path.ref'])
    queryj = pd.read_json(row['path.query'])

    ref_vars = pd.DataFrame(columns=['chrom','pos','ref','alt','variant'])
    query_vars = pd.DataFrame(columns=['chrom','pos','ref','alt','variant'])

    if 'columns' in refj['VARIANTS']['TIER1-3'].keys():
        ref_vars = pd.DataFrame(columns=refj['VARIANTS']['TIER1-3']['columns'],data=refj['VARIANTS']['TIER1-3']['data'])
        ref_vars['variant'] = ref_vars.apply(lambda x: x['gene']+':'+x['psyntax']+'['+x['vaf']+']',axis=1)

    if 'columns' in queryj['VARIANTS']['TIER1-3'].keys():
        query_vars = pd.DataFrame(columns=queryj['VARIANTS']['TIER1-3']['columns'],data=queryj['VARIANTS']['TIER1-3']['data'])
        query_vars['variant'] = query_vars.apply(lambda x: x['gene']+':'+x['psyntax']+'['+x['vaf']+']',axis=1)


    var_union = pd.merge(ref_vars[['chrom','pos','ref','alt','variant']],
                         query_vars[['chrom','pos','ref','alt','variant']],
                         on=['chrom','pos','ref','alt'],how='outer',suffixes=['.ref','.query'])

    cond = [
        var_union['variant.ref'].notna() & var_union['variant.query'].notna(),
        var_union['variant.ref'].isna() & var_union['variant.query'].notna(),
        var_union['variant.ref'].notna() & var_union['variant.query'].isna()
    ]

    vals = [
        'TP',
        'FP',
        'FN'
    ]

    var_union['type'] = pd.Categorical(np.select(cond, vals, default='Unknown'),categories=['TP','FP','FN'])

    out = dict(zip(out.columns.tolist(),[row['case'],
                                         os.path.abspath(row['path.ref']),
                                         os.path.abspath(row['path.query']),
                                         ref_vars.shape[0],query_vars.shape[0],
                                         ';'.join(var_union[var_union['type']=='FN']['variant.ref']),
                                         ';'.join(var_union[var_union['type']=='FP']['variant.query'])]))
    out['TP'] = var_union['type'].value_counts().to_dict()['TP']
    out['FP'] = var_union['type'].value_counts().to_dict()['FP']
    out['FN'] = var_union['type'].value_counts().to_dict()['FN']
    
    if out['reference_reportable_variants'] > 0:
        out['recall'] = round(out['TP']/(out['TP']+out['FN'])*100,1)
    else:
        out['recall'] = pd.NA

    if out['TP']+out['FP'] > 0:
        out['precision'] = round(out['TP']/(out['TP']+out['FP'])*100,1)
    else:
        out['precision'] = pd.NA

    if  out['precision'] is not pd.NA and out['recall'] is not pd.NA:
        out['f1'] = round(2 * (out['precision'] * out['recall']) / (out['precision'] + out['recall']),1)
    else:
        out['f1'] = pd.NA

    return(out)

#
# Script
#

parser = argparse.ArgumentParser(description='Compare MyeloSeqHD cases')
parser.add_argument('referenceDir',help='Reference case dir')
parser.add_argument('queryDir',help='Query case dir')
parser.add_argument('-o','--out',default=None,help='Query case dir')

args = parser.parse_args()

outfile = sys.stdout
if args.out is not None:
    outfile = args.out

refDat = pd.DataFrame()
refDat['path'] = search_files(args.referenceDir, "*.report.json")
refDat['case'] = refDat.apply(lambda x: os.path.splitext(os.path.basename(x['path']))[0].replace(".report",""),axis=1)

queryDat = pd.DataFrame()
queryDat['path'] = search_files(args.queryDir, "*.report.json")
queryDat['case'] = queryDat.apply(lambda x: os.path.splitext(os.path.basename(x['path']))[0].replace(".report",""),axis=1)

dat = pd.merge(refDat,queryDat,on='case',suffixes=['.ref','.query'],how='inner')

if dat.shape[0] < 1:
    sys.exit("No reports for matching cases in reference and query directories")


compared = dat.apply(lambda row: compare_cases(row),axis=1,result_type='expand')[['case','reference_reportable_variants','query_reportable_variants','TP','FN','FP','precision','recall','f1','reference_unique_reportable_variants','query_unique_reportable_variants','reference','query']]

compared.to_csv(outfile,sep="\t",index=False,na_rep='NA')

