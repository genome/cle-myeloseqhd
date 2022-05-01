#!/usr/bin/env python3

import os
import re
import sys
import csv
import glob
import time
import sqlite3
import fnmatch
import argparse
import warnings

from datetime import datetime
from cyvcf2 import VCF, Writer

def parse_csq_header(vcf_file):
    for header in vcf_file.header_iter():
        info = header.info(extra=True)
        if b'ID' in info.keys() and info[b'ID'] == b'CSQ':
            format_pattern = re.compile('Format: (.*)"')
            match = format_pattern.search(info[b'Description'].decode())
            return match.group(1).split('|')

def get_transcript_list(bed):
    t_list = []
    with open(bed, newline = '') as qcbed:
        bed_reader = csv.reader(qcbed, delimiter='\t')
        for line in bed_reader:
            t_list.append(line[8])
    qcbed.close()
    setlist = set(t_list) 
    return list(setlist)


mrn, acn, in_mrn, in_acn, db, vcf, path, covqc_bed = '', '', '', '', '', '', '', ''

parser = argparse.ArgumentParser(description='Upload to MyeloseqHD variantDB and query it')
parser.add_argument('-a',"--accession_query",type=str,help='Query DB with sample accession id')
parser.add_argument('-c',"--covqcbed",type=str,help='CoverageQC bed file for trusted gene/transcript to pick from vep annotation')
parser.add_argument('-d',"--db",type=str,help='Path to MyeloseqHD varaint sqlite3 DB file')
parser.add_argument('-i',"--mrn_input",type=str,help='Input MRN number for DB upload')
parser.add_argument('-j',"--accession_input",type=str,help='Input accession id for DB upload')
parser.add_argument('-m',"--mrn_query",type=str,help='Query DB with sample MRN number')
parser.add_argument('-p',"--path",type=str,help='Path to vcf files')
parser.add_argument('-v',"--vcf",type=str,help='Path to vcf file with variants')

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

if args.accession_query:
    acn = args.accession_query
if args.covqcbed:
    covqc_bed = args.covqcbed
if args.db:
    db = args.db
if args.mrn_input:
    in_mrn = args.mrn_input
if args.accession_input:
    in_acn = args.accession_input
if args.mrn_query:
    mrn = args.mrn_query
if args.path:
    path = args.path
if args.vcf:
    vcf = args.vcf

if not db:
    sys.exit("Provide variant sqlite3 database file")
if path and vcf:
    sys.exit("Provide either vcf directory path or vcf file")
if not os.path.exists(db) and mrn and acn:
    sys.exit("Need a valid varaint sqlite  DB to make queries")
if vcf or path:        
    if not os.path.exists(covqc_bed) or os.path.getsize(covqc_bed) == 0:
        sys.exit("Need a valid myeloseq CoverageQC bed file for helping vep annotation pick")
if vcf:
    if not os.path.exists(vcf) or os.path.getsize(vcf) == 0:
        sys.exit("Need a valid VCF file to upload variants")
if in_mrn and in_acn and path:
    sys.exit("Input MRN and input accession are for single vcf only, not vcf path")
if in_mrn and in_acn and vcf:
    if in_mrn == 'NONE' or in_acn == 'NONE':
        quit()

con = sqlite3.connect(db)
cur = con.cursor()

if vcf or path:   #upload db
    if os.path.exists(db) and os.path.getsize(db) > 0:
        print("MyeloseqHD varaint DB existing")
    else:
        print("Create MyeloseqHD variant table")
        cur.execute("create table myeloseqhd_variants (id integer primary key, mrn text, accession text, date text, version text, chromosome text, position integer, reference text, variant text, filter text, transcript_name text, consequence text, symbol text, gene_id text, exon text, intron text, p_syntax text, c_syntax text, coverage integer, vaf real, tamp integer, samp integer, amps text, UNIQUE (mrn, accession, chromosome, position, reference, variant) ON CONFLICT IGNORE)")
   
    if fnmatch.fnmatch(os.path.basename(covqc_bed), '*.b37.*'):
        version = 'v1'
    elif fnmatch.fnmatch(os.path.basename(covqc_bed), '*.hg38.*'):
        version = 'v2'
    else:
        sys.exit("Unrecognized coverage QC bed: " + covqc_bed)

    if vcf:
        all_vcf = [vcf]
    elif path:
        all_vcf = sorted(glob.glob(path + "/*.annotated_filtered.vcf.gz", recursive=True))
    else:
        sys.ext("provide either vcf directory or vcf file path")

    for f in all_vcf:
        if in_mrn and in_acn:
            new_mrn, new_acn = in_mrn, in_acn
        else:
            pattern = r"[A-Z]{4}\-(\d+)\-([A-Z]\d+\-\d+)"
            result = re.search(pattern, os.path.basename(f))
            if result:
                new_mrn, new_acn = result.group(1), result.group(2)
            else:
                warnings.warn("SKIP. Fail to get mrn and accession for: " + f)
                continue

        if version == 'v1':
            pattern = r"_(\d{4}\-\d{2}\-\d{2}_\d{2}:\d{2}:\d{2})\."
            result = re.search(pattern, os.path.basename(f))
            if result:
                date = result.group(1)
            else:
                sys.exit("Fail to get date for: " + f)
        elif version == 'v2':
            date = time.strftime('%Y-%m-%d_%H:%M:%S', time.localtime(os.path.getmtime(f)))

        vcf_parser = VCF(f)
        csq_fields = parse_csq_header(vcf_parser)
        trs_list   = get_transcript_list(covqc_bed)

        for variant in vcf_parser:
            #cyvcf variant.FILTER is None if variant FILTER column is PASS
            if variant.FILTER and not variant.INFO.get('MyeloSeqHDDB'):
                continue

            filter_str = 'PASS'
            if variant.FILTER is not None:
                filter_str = str(variant.FILTER).replace(";", "_") #semicolon and comma bring trouble to INFO tag

            chrom, pos, ref, var = str(variant.CHROM), variant.POS, str(variant.REF), str(variant.ALT[0])
            tamp, samp, amps     = int(variant.format("TAMP")[0][0]), int(variant.format("SAMP")[0][0]), variant.format("AMPS")[0]
            if not amps:
                amps = 'N/A'

            if version == 'v1':
                cov = int(variant.format("NR")[0][0])
                vaf = round(float(variant.format("CVAF")[0][0]), 3)
            elif version == 'v2':
                cov = int(variant.format("DP")[0][0])
                vaf = round(float(variant.format("VAF")[0][0]), 3)

            csq = variant.INFO.get('CSQ')
            if csq is None:
                sys.exit("No VEP fields: " + f)
            
            valid_annots = []
            pick_annot = {}

            for entry in csq.split(','):
                values = entry.split('|')
                annot = {}
                for key, value in zip(csq_fields, values):
                    annot[key] = value
                if annot['PICK'] == '1':
                    if annot['Feature'] in trs_list:
                        pick_annot = annot.copy()
                        break
                else:
                    if annot['Feature'] in trs_list:
                        valid_annots.append(annot.copy())

            if not pick_annot:
                if valid_annots:
                    pick_annot = valid_annots[0]
                else:
                    sys.exit("Failed to find valid csq annotation")

            gene_id, gene_name, transcript = pick_annot['Gene'], pick_annot['SYMBOL'], pick_annot['Feature']
            consequence, exon, intron      = pick_annot['Consequence'], pick_annot['EXON'], pick_annot['INTRON']
            
            if not exon:
                exon = 'N/A'
            if not intron:
                intron = 'N/A'

            csyntax = pick_annot['HGVSc'].split(":")
            if len(csyntax) > 1:
                c_syntax = csyntax[1]
            else:
                c_syntax = 'noncoding'

            psyntax = pick_annot['HGVSp'].split(":")
            if len(psyntax) > 1:
                p_syntax = psyntax[1]
            else:
                p_syntax = 'splice'

            cur.execute("insert into myeloseqhd_variants (mrn, accession, date, version, chromosome, position, reference, variant, filter, transcript_name, consequence, symbol, gene_id, exon, intron, p_syntax, c_syntax, coverage, vaf, tamp, samp, amps) values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (new_mrn, new_acn, date, version, chrom, pos, ref, var, filter_str, transcript, consequence, gene_name, gene_id, exon, intron,  p_syntax, c_syntax, cov, vaf, tamp, samp, amps))
            con.commit()

elif mrn and acn:   #query db and output a vcf file
    uniq_acn = cur.execute("select distinct accession from myeloseqhd_variants where mrn = ? and accession != ?", (mrn, acn)).fetchall()
    acn_list = []
    for row_acn in uniq_acn:
        acn_list.append(row_acn[0])

    var_rows = cur.execute("select * from myeloseqhd_variants where mrn = ? and accession != ?", (mrn, acn)).fetchall()
    desc_list = [i[0] for i in cur.description]
    desc_list = [e for e in desc_list if e not in ('id', 'chromosome', 'position', 'reference', 'variant', 'amps')]

    db_version = 'v1'
    outvcf_file = mrn + '_' + acn + '_query.vcf'
    vcf_writer = Writer.from_string(outvcf_file, "##fileformat=VCFv4.2\n#" + "\t".join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','MYELOSEQHDDB']),mode="w")
    vcf_writer.add_to_header('##MyeloSeqHDDB_VERSION=' + db_version)
    
    if acn_list:
        vcf_writer.add_to_header('##MyeloSeqHDDB_ACCESSIONLIST=' + ",".join(acn_list))
    
    vcf_writer.add_to_header('##MyeloSeqHDDB_FIELDS=' + "|".join(desc_list))
    vcf_writer.add_to_header('##contig=<ID=chr1,length=248956422>')
    vcf_writer.add_to_header('##contig=<ID=chr2,length=242193529>')
    vcf_writer.add_to_header('##contig=<ID=chr3,length=198295559>')
    vcf_writer.add_to_header('##contig=<ID=chr4,length=190214555>')
    vcf_writer.add_to_header('##contig=<ID=chr5,length=181538259>')
    vcf_writer.add_to_header('##contig=<ID=chr6,length=170805979>')
    vcf_writer.add_to_header('##contig=<ID=chr7,length=159345973>')
    vcf_writer.add_to_header('##contig=<ID=chr8,length=145138636>')
    vcf_writer.add_to_header('##contig=<ID=chr9,length=138394717>')
    vcf_writer.add_to_header('##contig=<ID=chr10,length=133797422>')
    vcf_writer.add_to_header('##contig=<ID=chr11,length=135086622>')
    vcf_writer.add_to_header('##contig=<ID=chr12,length=133275309>')
    vcf_writer.add_to_header('##contig=<ID=chr13,length=114364328>')
    vcf_writer.add_to_header('##contig=<ID=chr14,length=107043718>')
    vcf_writer.add_to_header('##contig=<ID=chr15,length=101991189>')
    vcf_writer.add_to_header('##contig=<ID=chr16,length=90338345>')
    vcf_writer.add_to_header('##contig=<ID=chr17,length=83257441>')
    vcf_writer.add_to_header('##contig=<ID=chr18,length=80373285>')
    vcf_writer.add_to_header('##contig=<ID=chr19,length=58617616>')
    vcf_writer.add_to_header('##contig=<ID=chr20,length=64444167>')
    vcf_writer.add_to_header('##contig=<ID=chr21,length=46709983>')
    vcf_writer.add_to_header('##contig=<ID=chr22,length=50818468>')
    vcf_writer.add_to_header('##contig=<ID=chrX,length=156040895>')
    vcf_writer.add_to_header('##contig=<ID=chrY,length=57227415>')

    vcf_writer.add_info_to_header({'ID': 'MyeloSeqHDDB','Description':'Comma separated list of variant info for the same MRN from the MyeloSeqHDDB','Type': 'String', 'Number': '.'})
    vcf_writer.add_format_to_header({'ID': 'GT','Description':'Genotype','Type': 'String', 'Number': '1'})
    vcf_writer.write_header()

    for r in var_rows:
        #id, mrn, accession, date, version, chromosome, position, reference, variant, filter, transcript_name, consequence, symbol, gene_id, exon, intron, p_syntax, c_syntax, coverage, vaf, tamp, samp, amps
        info = "|".join([r[1], r[2], r[3], r[4], r[9], r[10], r[11], r[12], r[13], r[14], r[15], r[16], r[17], str(r[18]), str(r[19]), str(r[20]), str(r[21])])
        new_rec = vcf_writer.variant_from_string("\t".join([r[5], str(r[6]), '.', r[7], r[8], '.', 'PASS', 'MyeloSeqHDDB='+info, 'GT', '0/1']))
        vcf_writer.write_record(new_rec)
    vcf_writer.close()

    #csv_path = mrn + '_' + acn + '_' + "query.tsv"
    #with open(csv_path, "w") as csv_fh:
    #    csv_writer = csv.writer(csv_fh, delimiter = "\t")
    #    csv_writer.writerow([i[0] for i in cur.description])
    #    csv_writer.writerows(cur)
    #csv_fh.close()
else:
    sys.exit("Either upload or query variant DB")

con.close()
