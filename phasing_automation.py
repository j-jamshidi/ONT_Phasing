import pysam
import argparse
import subprocess
from collections import Counter
import whatshap
from whatshap.vcf import VcfReader, VariantTable
from whatshap.align import edit_distance
from whatshap.variants import ReadSetReader
import os
import sys
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
    
    # Handle SNVs 
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

    # Run WhatsHap phase with output capture
    phase_command = f"whatshap phase -o {phased_vcf} --reference {reference_fasta} --ignore-read-groups {vcf_file} {bam_file}"
    print("Running WhatsHap phase command:")
    print(phase_command)
    
    try:
        phase_result = subprocess.run(phase_command, shell=True, check=True, capture_output=True, text=True)
        print("\nWhatsHap Phase stdout:")
        print(phase_result.stdout)
        
        if phase_result.stderr:
            print("\nWhatsHap Phase stderr:")
            print(phase_result.stderr)
    except subprocess.CalledProcessError as e:
        print("\nError running WhatsHap phase:")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise

    # Compress the phased VCF file with BGZIP
    bgzip_command = f"bgzip -f {phased_vcf}"
    subprocess.run(bgzip_command, shell=True, check=True)

    # Index the compressed phased VCF file
    index_command = f"tabix -p vcf {phased_vcf_gz}"
    subprocess.run(index_command, shell=True, check=True)

    # Tag reads with phase information
    haplotag_command = f"whatshap haplotag --tag-supplementary -o {output_bam} --reference {reference_fasta} {phased_vcf_gz} {bam_file} --ignore-read-groups"
    print("\nRunning WhatsHap haplotag command:")
    print(haplotag_command)
    
    try:
        haplotag_result = subprocess.run(haplotag_command, shell=True, check=True, capture_output=True, text=True)
        print("\nWhatsHap Haplotag stdout:")
        print(haplotag_result.stdout)
        
        if haplotag_result.stderr:
            print("\nWhatsHap Haplotag stderr:")
            print(haplotag_result.stderr)
    except subprocess.CalledProcessError as e:
        print("\nError running WhatsHap haplotag:")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        raise

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

    print("\nQuality Control:")

    for read in input_bam.fetch():
        if read.reference_start <= var1.pos - 1 and read.reference_end >= var2.pos:
            spanning_reads += 1

            alleles = [get_allele(read, variant) for variant in variants]

            if alleles[0] in ['ref', 'alt'] and alleles[1] in ['ref', 'alt']:
                clean_spanning_reads += 1

                if read.mapping_quality >= 15:
                    high_quality_reads += 1
                    output_bam.write(read)

    input_bam.close()
    output_bam.close()

    print("\nVariants to be phased:")
    print(f"Variant one: {var1.chrom} {var1.pos} {var1.ref} {var1.alts[0]}")
    print(f"Variant two: {var2.chrom} {var2.pos} {var2.ref} {var2.alts[0]}")
    print(f"The distance between variants is {var2.pos - var1.pos} bp")

    print(f"\nTotal reads: {total_reads}")
    print(f"Spanning reads: {spanning_reads} ({spanning_reads/total_reads*100:.2f}%)")
    print(f"Clean spanning reads: {clean_spanning_reads} ({clean_spanning_reads/spanning_reads*100:.2f}% of spanning reads)")
    print(f"High quality spanning reads (MAPQ >= 15): {high_quality_reads} ({high_quality_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)")

    if high_quality_reads > 50:
        print("QC PASSED for the quality (Q>15) and depth (>50) of sequencing")
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

    for read in bam.fetch():
        alleles = [get_allele(read, variant) for variant in variants]

        if alleles[0] == 'ref' and alleles[1] == 'ref':
            ref_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'alt':
            alt_reads += 1
        elif alleles[0] == 'ref' and alleles[1] == 'alt':
            ref_alt_reads += 1
        elif alleles[0] == 'alt' and alleles[1] == 'ref':
            alt_ref_reads += 1

    cis_reads = ref_reads + alt_reads
    trans_reads = ref_alt_reads + alt_ref_reads
    cis_percentage = cis_reads / total_reads * 100 if total_reads > 0 else 0
    trans_percentage = trans_reads / total_reads * 100 if total_reads > 0 else 0
    print("Result:\n")
    print(f"\nTotal high quality spanning reads: {total_reads}")
    print("\nDetailed categorisation:")
    print(f"Reads with ref allele for both variants (Cis): {ref_reads} ({ref_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt allele for both variants (Cis): {alt_reads} ({alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with ref for first, alt for second (Trans): {ref_alt_reads} ({ref_alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt for first, ref for second (Trans): {alt_ref_reads} ({alt_ref_reads/total_reads*100:.2f}%)")
    print(f"\nCIS reads (both ref or both alt): {cis_reads} ({cis_percentage:.2f}%)")
    print(f"TRANS reads (one ref, one alt): {trans_reads} ({trans_percentage:.2f}%)")

    if cis_percentage > trans_percentage:
        print(f"\nCounting the reads determined the phase as CIS \n\nChimeric reads percentage: {trans_percentage:.2f}% ")
    else:
        print(f"Counting reads determined the phase as TRANS \n\nChimeric reads percentage: {cis_percentage:.2f}% ")

    bam.close()

def parse_vcf_phasing(phased_vcf_file):
    """
    Read the phased VCF file and determine if variants are in CIS or TRANS.
    
    Args:
        phased_vcf_file (str): Path to the phased VCF file
    
    Returns:
        str: Interpretation of variant phasing (CIS or TRANS)
    """
    vcf = pysam.VariantFile(phased_vcf_file)
    variants = list(vcf.fetch())
    
    if len(variants) != 2:
        raise ValueError("VCF file should contain exactly two variants")
    
    # Extract genotype information
    var1_genotype = variants[0].samples[0]['GT']
    var2_genotype = variants[1].samples[0]['GT']
    
    # Determine phasing
    if var1_genotype == var2_genotype:
        return "CIS"
    else:
        return "TRANS"

def main():
    """
    Main function to run the complete analysis pipeline.
    """
    # Parse command-line arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.bam> <input.vcf>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]

    # Create output names
    output_bam = "phased_" + os.path.basename(bam_file)
    output_txt = "variant_phasing_report.txt"

    # Redirect stdout and stderr to the output text file
    with open(output_txt, 'w') as output_file:
        # Redirect standard output to the file
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = output_file
        sys.stderr = output_file

        try:
            # Quality Control
            clean_span_bam = "clean-span-hq.bam"
            high_quality_reads = quality_control(bam_file, vcf_file, clean_span_bam)

            # Index the output BAM file
            index_bam_command = f"samtools index {clean_span_bam}"
            subprocess.run(index_bam_command, shell=True, check=True)

            if high_quality_reads > 0:
                print("\nProceeding with phasing and analysis...")
                
                # Capture WhatsHap output
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

                # Determine phasing from VCF file
                phased_vcf = os.path.splitext(phased_vcf_gz)[0]
                vcf_phasing = parse_vcf_phasing(phased_vcf)
                
                
            
                print(" "*40)
                print(f"WhatsHap analysis determined the phase as {vcf_phasing}")
            else:
                print("No high quality spanning reads found. Analysis cannot proceed.")

        except Exception as e:
            print(f"An error occurred: {e}")
        finally:
            # Restore stdout and stderr
            sys.stdout = original_stdout
            sys.stderr = original_stderr

    print(f"Analysis complete. Results written to {output_txt}")

if __name__ == "__main__":
    main()