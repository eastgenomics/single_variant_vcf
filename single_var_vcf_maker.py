#!/usr/bin/env python3

import re
import sys
import pandas as pd

df = pd.read_excel(sys.argv[1], dtype=str, engine="odf")

# Define static elements of the VCF
B37_HEADER = """##fileformat=VCFv4.2
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

B38_HEADER = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=V1_evidence,Number=.,Type=String,Description="Free text from scientists">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##contig=<ID=chr1,length=248956422,assembly=hg38>
##contig=<ID=chr2,length=242193529,assembly=hg38>
##contig=<ID=chr3,length=198295559,assembly=hg38>
##contig=<ID=chr4,length=190214555,assembly=hg38>
##contig=<ID=chr5,length=181538259,assembly=hg38>
##contig=<ID=chr6,length=170805979,assembly=hg38>
##contig=<ID=chr7,length=159345973,assembly=hg38>
##contig=<ID=chr8,length=145138636,assembly=hg38>
##contig=<ID=chr9,length=138394717,assembly=hg38>
##contig=<ID=chr10,length=133797422,assembly=hg38>
##contig=<ID=chr11,length=135086622,assembly=hg38>
##contig=<ID=chr12,length=133275309,assembly=hg38>
##contig=<ID=chr13,length=114364328,assembly=hg38>
##contig=<ID=chr14,length=107043718,assembly=hg38>
##contig=<ID=chr15,length=101991189,assembly=hg38>
##contig=<ID=chr16,length=90338345,assembly=hg38>
##contig=<ID=chr17,length=83257441,assembly=hg38>
##contig=<ID=chr18,length=80373285,assembly=hg38>
##contig=<ID=chr19,length=58617616,assembly=hg38>
##contig=<ID=chr20,length=64444167,assembly=hg38>
##contig=<ID=chr21,length=46709983,assembly=hg38>
##contig=<ID=chr22,length=50818468,assembly=hg38>
##contig=<ID=chrX,length=156040895,assembly=hg38>
##contig=<ID=chrY,length=57227415,assembly=hg38>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"""

HEADER_DICT = {"GRCh37":B37_HEADER, "GRCh38":B38_HEADER}

ID = QUAL = FILTER = "."
FORMAT = "GT:AD"
AD = "0,0"

# Regex pattern for extracting chrom,pos,ref,alt
PATTERN = r"([1-9]?|1[0-9]?|2[0-2]?|X|Y):(\d+):([ACGT]+):([ACGT]+)"


# Define functions

def get_variant(r_index, row_data):
    """Parses and checks variant information from input"""
    result = re.search(PATTERN, row_data)
    if result is not None:
        if len(result[3]) > 1 and len(result[4]) > 1:
            raise ValueError(
                f"REF/ALT are not valid for row index: {r_index}") 
        return result
    raise ValueError(
        f"Variant notation is not valid for row index: {r_index}")

def get_zygosity(r_index, row_data):
    """Parses and checks zygosity from input"""
    if row_data in ["het", "Het", "heterozygous", "Heterozygous"]:
        return "0/1"
    if row_data in ["hom", "Hom", "homozygous", "Homozygous"]:
        return "1/1"
    raise ValueError(f"Zygosity is not valid for row index: {r_index}")
    
# Loop and parse out variants in the excel and export in VCF format
for index, row in df.iterrows():
    
    g_build = row["Genome_build"].strip()
    # At the moment only accepting b37
    if g_build in ["GRCh37"]:    
        vcf_header = HEADER_DICT[g_build]
    else:
        raise ValueError(
            f"The genome build is not valid for row index: {index}")

    if pd.isna(row["Participant_ID"]):
        raise ValueError(
            f"The sample ID could not be found for row index: {index}")

    sample_id = row["Participant_ID"].strip()
    vcf_header = f"{vcf_header}\t{sample_id}\n"

    variant = get_variant(index, row["V1"].strip())
    chrom = variant[1]
    pos = variant[2]
    ref = variant[3]
    alt = variant[4]
    
    GT = get_zygosity(index, row["V1_zygosity"].strip())
    gt_ad = f"{GT}:{AD}"

    info = 'V1_evidence="{}"'.format(row["V1_evidence"])

    RECORD = "\t".join([chrom, pos, ID, ref, alt, QUAL, FILTER, info,
                        FORMAT, gt_ad])

    output_vcf = f"{sample_id}.vcf"
    with open(output_vcf, "w", encoding="UTF-8") as vcf:
        vcf.write(vcf_header+RECORD)
        