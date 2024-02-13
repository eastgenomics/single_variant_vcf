"""Generate VCF for each variant in spreadsheet"""
import re
import argparse
import pandas as pd

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

HEADER_DICT = {"GRCh37": B37_HEADER, "GRCh38": B38_HEADER}

ID = QUAL = FILTER = "."
FORMAT = "GT:AD"
AD = "0,0"

# Regex pattern for extracting chrom,pos,ref,alt
PATTERN = r"([1-9]?|1[0-9]?|2[0-2]?|X|Y):(\d+):([ACGT]+):([ACGT]+)"


def parse_args():
    """
    Parse arguments given at cmd line.

    Returns:
        - args (Namespace): object containing parsed arguments.
    """

    parser = argparse.ArgumentParser(
        description="Generate VCF file for each variant"
    )

    parser.add_argument(
        '-i', '--input',
        help="Name of input spreadsheet containing variant information",
        required=True,
        type=str
    )

    parser.add_argument(
        '-g', '--genome',
        help="Name of column in the input spreadsheet containing genome build",
        default="Genome_build",
        type=str
    )

    parser.add_argument(
        '-p', '--participant',
        help="Name of column in the input spreadsheet containing participant\
            ID",
        default="Participant_ID",
        type=str
    )

    parser.add_argument(
        '-v', '--variant',
        help="Name of column in the input spreadsheet containing variant\
            description",
        default="V1",
        type=str
    )

    parser.add_argument(
        '-z', '--zygosity',
        help="Name of column in the input spreadsheet containing zygosity\
            description",
        default="V1_zygosity",
        type=str
    )

    parser.add_argument(
        '-e', '--evidence',
        help="Name of column in the input spreadsheet containing\
            additional information for the assertion of pathogenicity",
        default="V1_evidence",
        type=str
    )

    args = parser.parse_args()
    return args


def get_variant(r_index: int, row_data: pd.Series) -> re.Match:
    """
    Parses and checks variant information from input

    Raises:
        ValueError: if REF/ALT values are invalid
        ValueError: if variant notation is invalid

    Returns:
        re.Match: Regex match object containing chrom, pos, ref and alt values
    """
    result = re.search(PATTERN, row_data)
    if result is not None:
        if len(result[3]) > 1 and len(result[4]) > 1:
            raise ValueError(
                f"REF/ALT are not valid for row index: {r_index}")
        return result
    raise ValueError(
        f"Variant notation is not valid for row index: {r_index}")


def get_zygosity(r_index: int, row_data: pd.Series) -> str:
    """
    Parses and checks zygosity from input

    Raises:
        ValueError: if zygosity is note valid

    Returns:
        str: zygosity value
    """
    if row_data.lower() in ["heterozygous", "het"]:
        return "0/1"
    if row_data.lower() in ["homozygous", "hom"]:
        return "1/1"
    raise ValueError(f"Zygosity is not valid for row index: {r_index}")


def make_vcf(df: pd.DataFrame, genome_build_col: str, participant_id_col: str,
             variant_col: str, zygosity_col: str, evidence_col: str):
    """
    Loop over and parse out variants in the excel and export in VCF format
    Args:
        df (pd.DataFrame): dataframe object of input spreadsheet
        participant_id_col (str): name of column containing participant IDs
        variant_col (str): name of column containing variant description
        zygosity_col (str): name of column containing zygosity description
        evidence_col (str): name of column containing additional information
            for the assertion of pathogenicity

    Raises:
        ValueError: if genome build is not valid
        ValueError: if participant ID cannot not be found

    Outputs:
        VCF file containing single variant
    """
    for index, row in df.iterrows():
        genome_build = row[genome_build_col].strip()
        if genome_build in HEADER_DICT.keys():
            vcf_header = HEADER_DICT[genome_build]
        else:
            raise ValueError(
                f"The genome build is not valid for row index: {index}"
            )

        if pd.isna(row[participant_id_col]):
            raise ValueError(
                f"The participant ID could not be found for row index: {index}"
            )

        participant_id = row[participant_id_col].strip()
        vcf_header = f"{vcf_header}\t{participant_id}\n"

        variant = get_variant(index, row[variant_col].strip())
        chrom = variant[1]
        pos = variant[2]
        ref = variant[3]
        alt = variant[4]

        gt = get_zygosity(index, row[zygosity_col].strip())
        gt_ad = f"{gt}:{AD}"

        info = f'V1_evidence="{row[evidence_col]}"'

        record = "\t".join([chrom, pos, ID, ref, alt, QUAL, FILTER, info,
                            FORMAT, gt_ad])

        output_vcf = f"{participant_id}_{genome_build}.vcf"
        with open(output_vcf, "w", encoding="UTF-8") as vcf:
            vcf.write(vcf_header+record)


def main():
    """Main function to generate VCF for each variant in spreadsheet"""
    args = parse_args()
    df = pd.read_excel(args.input, dtype=str, engine="odf")
    make_vcf(
        df=df,
        genome_build_col=args.genome,
        participant_id_col=args.participant,
        variant_col=args.variant,
        zygosity_col=args.zygosity,
        evidence_col=args.evidence
    )


if __name__ == "__main__":
    main()
