#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import re
import sys

# Read in input excel
df = pd.read_excel(sys.argv[1], dtype=str, engine='odf')

# Define static elements of the VCF
static_header = '''##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=V1_evidence,Number=.,Type=String,Description="Free text from scientists">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT'''
id_val = qual = filter_val = "."
format_val = "GT:AD"
ad = "0,0"

# Regex pattern for extracting chrom,pos,ref,alt  
var_pattern = "(1[0-9]?|2[0-2]?|X|Y):(\d+):([ACGT]+):([ACGT]+)"

# Loop and parse out variants in the excel and export in VCF format
for index,row in df.iterrows():
    
    if pd.isna(row['Participant_ID']):
        print("The sample ID could not be found for row index: {}".format(index))
        continue
    
    sample_id = row['Participant_ID'].strip()
    vcf_header = "{}\t{}\n".format(static_header,sample_id)
    
    variant = row['V1'].strip()
    result = re.search(var_pattern,variant)
    try:
        chrom = result[1]
    except TypeError:
        print("CHROM could not be found for row index: {}".format(index))
        continue
    
    try:
        pos = result[2]
    except TypeError:
        print("POS could not be found for row index: {}".format(index))
        continue
    
    try:
        ref = result[3]
    except TypeError:
        print("REF could not be found for row index: {}".format(index))
        continue
    
    try:
        alt = result[4]
    except TypeError:
        print("ALT could not be found for row index: {}".format(index))
        continue
    
    if len(ref) > 1 and len(alt) > 1:
        print("REF/ALT are not valid for row index: {}".format(index))
        continue    
    
    zygos = row['V1_zygosity'].strip()
    if zygos in ["het", "Het", "heterozygous", "Heterozygous"]:
        gt = "0/1"
    elif zygos in ["hom", "Hom", "homozygous", "Homozygous"]:
        gt = "1/1"
    else:
        print("Zygosity is not valid for row index: {}".format(index))
        continue
    
    gt_ad = "{}:{}".format(gt,ad)
    
    info = "V1_evidence=\"{}\"".format(row['V1_evidence'])
    
    record = "\t".join([chrom,pos,id_val,ref,alt,qual,filter_val,info,format_val,gt_ad])
    
    output_VCF = "{}.vcf".format(sample_id)
    with open(output_VCF, 'w') as vcf:
        vcf.write(vcf_header+record)

    
    
    
        
        
        
        
    
    
  
    



