import pysam
import argparse
import subprocess
from collections import Counter
from tqdm import tqdm
import whatshap
from whatshap.vcf import VcfReader, VariantTable
from whatshap.align import edit_distance
import os

REFERENCE_FASTA = "/Users/javadjamshidi/Desktop/Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

def run_whatshap(bam_file, vcf_file, output_bam, reference_fasta):
    # Run WhatsHap phase
    phase_command = f"whatshap phase -o phased.vcf --reference {reference_fasta} --ignore-read-groups {vcf_file} {bam_file}"
    subprocess.run(phase_command, shell=True, check=True)

    # Compress the phased VCF file with BGZIP
    bgzip_command = "bgzip -f phased.vcf"
    subprocess.run(bgzip_command, shell=True, check=True)

    # Index the compressed phased VCF file
    index_command = "tabix -p vcf phased.vcf.gz"
    subprocess.run(index_command, shell=True, check=True)

    # Tag reads with phase information
    haplotag_command = f"whatshap haplotag --tag-supplementary -o {output_bam} --reference {reference_fasta} phased.vcf.gz {bam_file} --ignore-read-groups"
    subprocess.run(haplotag_command, shell=True, check=True)

    # Index the output BAM file
    index_bam_command = f"samtools index {output_bam}"
    subprocess.run(index_bam_command, shell=True, check=True)

def get_allele(read, variant):
    pos = variant.pos - 1  # 0-based position
    if pos < read.reference_start or pos >= read.reference_end:
        return None

    aligned_pairs = read.get_aligned_pairs()
    query_pos = None
    for qpos, rpos in aligned_pairs:
        if rpos == pos:
            query_pos = qpos
            break

    if query_pos is None:
        return None

    ref = variant.ref
    alt = variant.alts[0] if variant.alts else None

    # Extract the relevant portion of the read sequence
    read_seq = read.query_sequence[query_pos:query_pos + max(len(ref), len(alt))]

    # Compare the read sequence to both ref and alt alleles
    ref_distance = edit_distance(read_seq, ref)
    alt_distance = edit_distance(read_seq, alt)

    if ref_distance < alt_distance:
        return 'ref'
    elif alt_distance < ref_distance:
        return 'alt'
    else:
        return None

def quality_control(input_bam, vcf_file, output_bam):
    vcf = pysam.VariantFile(vcf_file)
    variants = list(vcf.fetch())
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")

    var1, var2 = variants

    input_bam = pysam.AlignmentFile(input_bam, "rb")
    output_bam = pysam.AlignmentFile(output_bam, "wb", template=input_bam)

    total_reads = input_bam.count()
    spanning_reads = 0
    clean_spanning_reads = 0
    high_quality_reads = 0

    print("\nQuality Control:")
    pbar = tqdm(total=total_reads, desc="Processing reads")

    for read in input_bam.fetch():
        pbar.update(1)
        if read.reference_start <= var1.pos - 1 and read.reference_end >= var2.pos:
            spanning_reads += 1

            alleles = [get_allele(read, variant) for variant in variants]

            if alleles[0] in ['ref', 'alt'] and alleles[1] in ['ref', 'alt']:
                clean_spanning_reads += 1
                
                if read.mapping_quality >= 20:
                    high_quality_reads += 1
                    output_bam.write(read)

    pbar.close()
    input_bam.close()
    output_bam.close()

    print("\nVariants to be phased:")
    print(f"Variant one: {var1.pos} {var1.ref} {var1.alts[0]}")
    print(f"Variant two: {var2.pos} {var2.ref} {var2.alts[0]}")

    print(f"\nTotal reads: {total_reads}")
    print(f"Spanning reads: {spanning_reads} ({spanning_reads/total_reads*100:.2f}%)")
    print(f"Clean spanning reads: {clean_spanning_reads} ({clean_spanning_reads/spanning_reads*100:.2f}% of spanning reads)")
    print(f"High quality reads (MAPQ >= 20): {high_quality_reads} ({high_quality_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)")

    if high_quality_reads > 100:
        print("QC PASSED for the quality and depth of sequencing")
    else:
        print("WARNING! The number of passed reads are below QC. Proceed with caution!")

    return high_quality_reads

def analyze_reads(bam_file, vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    variants = list(vcf.fetch())
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")

    var1, var2 = variants

    bam = pysam.AlignmentFile(bam_file, "rb")

    total_reads = bam.count()
    ref_reads = 0
    alt_reads = 0
    ref_alt_reads = 0
    alt_ref_reads = 0

    pbar = tqdm(total=total_reads, desc="Analyzing reads")

    for read in bam.fetch():
        pbar.update(1)
        alleles = [get_allele(read, variant) for variant in variants]

        if alleles[0] == 'ref' and alleles[1] == 'ref':
            ref_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'alt':
            alt_reads += 1
        elif alleles[0] == 'ref' and alleles[1] == 'alt':
            ref_alt_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'ref':
            alt_ref_reads += 1

    pbar.close()

    cis_reads = ref_reads + alt_reads
    trans_reads = ref_alt_reads + alt_ref_reads
    cis_percentage = cis_reads / total_reads * 100
    trans_percentage = trans_reads / total_reads * 100

    print(f"\nTotal high quality spanning reads: {total_reads}")
    print("\nDetailed categorization:")
    print(f"Reads with reference allele for both variants: {ref_reads} ({ref_reads/total_reads*100:.2f}%)")
    print(f"Reads with alternate allele for both variants: {alt_reads} ({alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with ref for first, alt for second: {ref_alt_reads} ({ref_alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt for first, ref for second: {alt_ref_reads} ({alt_ref_reads/total_reads*100:.2f}%)")
    print(f"\nCIS reads (both ref or both alt): {cis_reads} ({cis_percentage:.2f}%)")
    print(f"TRANS reads (one ref, one alt): {trans_reads} ({trans_percentage:.2f}%)")

    print("\nResult:")
    if cis_percentage > trans_percentage:
        print(f"The variants are likely in CIS with {trans_percentage:.2f}% chimeric reads\nPlease confirm with cheking the phased vcf and bam files")
    else:
        print(f"The variants are likely in TRANS with {cis_percentage:.2f}% chimeric reads\nPlease confirm with cheking the phased vcf and bam files")

    bam.close()

def main():
    bam_file = input("Please enter the path to the BAM file: ").strip().strip("'\"")
    vcf_file = input("Please enter the path to the VCF file: ").strip().strip("'\"")

    output_bam = "phased_" + os.path.basename(bam_file)

    # Quality Control
    clean_span_bam = "clean-span-hq.bam"
    high_quality_reads = quality_control(bam_file, vcf_file, clean_span_bam)

    # Index the output BAM file
    index_bam_command = f"samtools index {clean_span_bam}"
    subprocess.run(index_bam_command, shell=True, check=True)

    if high_quality_reads > 0:
        print("\nProceeding with phasing and analysis...")
        # Run WhatsHap on clean spanning reads
        print("\n" + "="*40)
        print("WhatsHap output:")
        print("="*40)
        run_whatshap(clean_span_bam, vcf_file, output_bam, REFERENCE_FASTA)
        print("="*40 + "\n")

        # Analyze phased reads
        analyze_reads(output_bam, vcf_file)
    else:
        print("No high quality spanning reads found. Analysis cannot proceed.")

if __name__ == "__main__":
    main()