import pysam
import argparse
import subprocess
from collections import Counter
from tqdm import tqdm
import whatshap
from whatshap.vcf import VcfReader, VariantTable
from whatshap.align import edit_distance
from whatshap.variants import ReadSetReader
import os
from typing import Optional, Tuple, List

#Change the path to your own reference file
REFERENCE_FASTA = "/Users/javadjamshidi/Desktop/Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

def normalize_allele(allele: str) -> str:
    """Normalize allele to uppercase and replace empty alleles with <DEL>."""
    if allele == "":
        return "<DEL>"
    return allele.upper()

def get_allele(read, variant) -> Optional[str]:
    """
    Determine which allele is present in a read at a given variant position.
    
    Args:
        read: pysam.AlignedSegment object
        variant: pysam.VariantRecord object
    
    Returns:
        'ref', 'alt', or None if the allele cannot be determined
    """
    pos = variant.pos - 1  # Convert to 0-based position
    
    # Check if the variant position is covered by the read
    if pos < read.reference_start or pos >= read.reference_end:
        return None
    
    # Get the reference and alternate alleles
    ref_allele = variant.ref.upper()
    alt_allele = variant.alts[0].upper() if variant.alts else ""
    
    # Get the aligned pairs with sequence
    aligned_pairs = read.get_aligned_pairs(with_seq=True)
    
    # Extract the variant region with context
    context_size = max(len(ref_allele), len(alt_allele)) + 10  # Increased context size
    variant_region = []
    ref_pos_idx = None
    
    # Find the variant position and collect context
    for idx, (qpos, rpos, base) in enumerate(aligned_pairs):
        if rpos == pos:
            ref_pos_idx = idx
            # Collect preceding context
            start_idx = max(0, idx - context_size//2)
            variant_region.extend(aligned_pairs[start_idx:idx])
            break
    
    if ref_pos_idx is None:
        return None
    
    # Collect following context
    end_idx = min(len(aligned_pairs), ref_pos_idx + context_size//2)
    variant_region.extend(aligned_pairs[ref_pos_idx:end_idx])
    
    # Extract sequences for comparison
    read_seq = ""
    ref_seq = ""
    read_positions = []
    ref_positions = []
    
    for qpos, rpos, base in variant_region:
        if qpos is not None and rpos is not None:
            read_seq += read.query_sequence[qpos]
            ref_seq += base.upper() if base else 'N'
            read_positions.append(qpos)
            ref_positions.append(rpos)
        elif qpos is None and rpos is not None:
            # Deletion in read
            read_seq += '-'
            ref_seq += base.upper() if base else 'N'
            read_positions.append(None)
            ref_positions.append(rpos)
        elif qpos is not None and rpos is None:
            # Insertion in read
            read_seq += read.query_sequence[qpos]
            ref_seq += '-'
            read_positions.append(qpos)
            ref_positions.append(None)
    
    # Handle different variant types
    if len(ref_allele) != len(alt_allele):
        # Indel case
        variant_start_idx = ref_positions.index(pos) if pos in ref_positions else -1
        if variant_start_idx == -1:
            return None
            
        # For insertions 
        if len(alt_allele) > len(ref_allele):
            # Extract the sequence around the insertion point
            insertion_length = len(alt_allele) - len(ref_allele)
            sequence_window = 5  # Consider surrounding context
            
            # Get the read sequence at and after the variant position
            read_variant_seq = ""
            current_idx = variant_start_idx
            insertion_bases = 0
            context_bases = 0
            total_bases_needed = insertion_length + sequence_window
            
            while current_idx < len(read_positions) and len(read_variant_seq) < total_bases_needed:
                if read_positions[current_idx] is not None:
                    read_variant_seq += read_seq[current_idx]
                    
                    # Count insertion bases (those without corresponding reference position)
                    if ref_positions[current_idx] is None:
                        insertion_bases += 1
                    else:
                        context_bases += 1
                current_idx += 1
            
            # Check if we have enough bases for comparison
            if len(read_variant_seq) >= insertion_length:
                # Calculate alignment scores against both alleles
                alt_score = 0
                ref_score = 0
                
                # Compare with alt allele
                for i in range(min(len(read_variant_seq), len(alt_allele))):
                    if read_variant_seq[i] == alt_allele[i]:
                        alt_score += 1
                
                # Compare with ref allele
                for i in range(min(len(read_variant_seq), len(ref_allele))):
                    if read_variant_seq[i] == ref_allele[i]:
                        ref_score += 1
                
                # Add weight to the insertion detection
                if insertion_bases >= insertion_length * 0.8:  # Allow for some sequencing errors
                    alt_score += 2  # Bonus for having the right number of inserted bases
                
                # Calculate match percentages
                alt_match_percent = alt_score / len(alt_allele)
                ref_match_percent = ref_score / len(ref_allele)
                
                # Make the final call
                if alt_match_percent >= 0.8 and insertion_bases >= insertion_length * 0.8:
                    return 'alt'
                elif ref_match_percent >= 0.8 and insertion_bases == 0:
                    return 'ref'
                
        # For deletions 
        else:
            deletion_length = len(ref_allele) - len(alt_allele)
            expected_ref_positions = list(range(pos, pos + len(ref_allele)))
            
            # Count missing bases in the read
            missing_bases = 0
            for exp_pos in expected_ref_positions:
                if exp_pos not in ref_positions or read_positions[ref_positions.index(exp_pos)] is None:
                    missing_bases += 1
            
            # Calculate match scores for both alleles
            ref_match_score = len(ref_allele) - missing_bases
            alt_match_score = missing_bases
            
            # Add additional check for partial deletions
            if missing_bases >= deletion_length * 0.8:  # Allow some mismatches
                return 'alt'
            elif ref_match_score >= len(ref_allele) * 0.8:  # Allow some mismatches
                return 'ref'
    
    # Handle SNVs (unchanged)
    else:
        query_pos = aligned_pairs[ref_pos_idx][0]
        if query_pos is None:
            return None
        
        read_base = read.query_sequence[query_pos]
        if read_base == ref_allele:
            return 'ref'
        elif read_base == alt_allele:
            return 'alt'
    
    return None
    

def run_whatshap(bam_file, vcf_file, output_bam, reference_fasta):
    """
    Run WhatsHap for phasing variants and tagging reads.
    """
    # Generate the output VCF file name based on the input VCF file name
    vcf_prefix = os.path.splitext(os.path.basename(vcf_file))[0]
    phased_vcf = f"{vcf_prefix}_Phased.vcf"
    phased_vcf_gz = f"{phased_vcf}.gz"

    # Run WhatsHap phase
    phase_command = f"whatshap phase -o {phased_vcf} --reference {reference_fasta} --ignore-read-groups {vcf_file} {bam_file}"
    subprocess.run(phase_command, shell=True, check=True)

    # Compress the phased VCF file with BGZIP
    bgzip_command = f"bgzip -f {phased_vcf}"
    subprocess.run(bgzip_command, shell=True, check=True)

    # Index the compressed phased VCF file
    index_command = f"tabix -p vcf {phased_vcf_gz}"
    subprocess.run(index_command, shell=True, check=True)

    # Tag reads with phase information
    haplotag_command = f"whatshap haplotag --tag-supplementary -o {output_bam} --reference {reference_fasta} {phased_vcf_gz} {bam_file} --ignore-read-groups"
    subprocess.run(haplotag_command, shell=True, check=True)

    # Index the output BAM file
    index_bam_command = f"samtools index {output_bam}"
    subprocess.run(index_bam_command, shell=True, check=True)

    return phased_vcf_gz

def quality_control(input_bam, vcf_file, output_bam):
    """
    Perform quality control on input BAM file and filter reads based on quality metrics.
    """
    if not os.path.exists(input_bam + '.bai'):
        pysam.index(input_bam)
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
    low_qv_reads = 0

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
                    # Check QV for each base of the variants
                    qv_pass = True
                    if read.query_qualities is not None:
                        for variant in variants:
                            pos = variant.pos - 1  # 0-based position
                            aligned_pairs = read.get_aligned_pairs()
                            query_pos = None
                            for qpos, rpos in aligned_pairs:
                                if rpos == pos:
                                    query_pos = qpos
                                    break

                            if query_pos is not None and read.query_qualities[query_pos] < 10:
                                qv_pass = False
                                break

                    if qv_pass:
                        high_quality_reads += 1
                        output_bam.write(read)
                    else:
                        low_qv_reads += 1

    pbar.close()
    input_bam.close()
    output_bam.close()

    print("\nVariants to be phased:")
    print(f"Variant one: {var1.chrom} {var1.pos} {var1.ref} {var1.alts[0]}")
    print(f"Variant two: {var2.chrom} {var2.pos} {var2.ref} {var2.alts[0]}")

    print(f"\nTotal reads: {total_reads}")
    print(f"Spanning reads: {spanning_reads} ({spanning_reads/total_reads*100:.2f}%)")
    print(f"Clean spanning reads: {clean_spanning_reads} ({clean_spanning_reads/spanning_reads*100:.2f}% of spanning reads)")
    print(f"High quality reads (MAPQ >= 20 and QV >= 10): {high_quality_reads} ({high_quality_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)")
    print(f"Reads removed due to low QV: {low_qv_reads} ({low_qv_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)")

    if high_quality_reads > 100:
        print("QC PASSED for the quality (Q>20) and depth (>100) of sequencing")
    else:
        print("WARNING! The number of passed reads are below QC. Proceed with caution!")

    return high_quality_reads

def analyze_reads(bam_file, vcf_file):
    """
    Analyze phased reads to determine variant relationships (CIS/TRANS).
    """
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
    cis_percentage = cis_reads / total_reads * 100 if total_reads > 0 else 0
    trans_percentage = trans_reads / total_reads * 100 if total_reads > 0 else 0

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
        print(f"The variants are likely in CIS with {trans_percentage:.2f}% chimeric reads\nPlease confirm with checking the phased vcf and bam files")
    else:
        print(f"The variants are likely in TRANS with {cis_percentage:.2f}% chimeric reads\nPlease confirm with checking the phased vcf and bam files")

    bam.close()

def main():
    """
    Main function to run the complete analysis pipeline.
    """
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
        phased_vcf_gz = run_whatshap(clean_span_bam, vcf_file, output_bam, REFERENCE_FASTA)
        print("="*40 + "\n")

        # Analyze phased reads
        analyze_reads(output_bam, vcf_file)

        # Unzip the phased VCF file
        unzip_command = f"gunzip -f {phased_vcf_gz}"
        subprocess.run(unzip_command, shell=True, check=True)
    else:
        print("No high quality spanning reads found. Analysis cannot proceed.")

if __name__ == "__main__":
    main()
