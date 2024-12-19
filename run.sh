#!/bin/bash

# Remove header from CSV file
tail -n +2 sample_sheet.csv > sample_sheet_nh
printf '\r' >> sample_sheet_nh

# Read the CSV file line by line
while IFS=, read -r Batch Barcode Episode Family Initial Ampliconsize Coordinate Gene Variant1 Variant2 Distance EpisodeWES remainder
do 
    {
        echo $Barcode
        # Create dummy VCF file for this barcode/episode
        cp dummy.vcf bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/sample/${Episode}/g" bam_pass/${Barcode}/${Episode}.vcf
        
        # Convert : > and space to - in Variant1 and Variant2
        Variant1=$(echo ${Variant1} | tr ':> ' '-')
        Variant2=$(echo ${Variant2} | tr ':> ' '-')
        
        IFS='-' read -r -a array1 <<< ${Variant1}
        sed -i -e "s/chrA/${array1[0]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/POS1/${array1[1]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/REF1/${array1[2]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/ALT1/${array1[3]}/g" bam_pass/${Barcode}/${Episode}.vcf
        unset array1

        IFS='-' read -r -a array2 <<< ${Variant2}
        sed -i -e "s/chrB/${array2[0]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/POS2/${array2[1]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/REF2/${array2[2]}/g" bam_pass/${Barcode}/${Episode}.vcf
        sed -i -e "s/ALT2/${array2[3]}/g" bam_pass/${Barcode}/${Episode}.vcf
        unset array2

        # Sort the VCF file after replacements
        grep -v "^#" bam_pass/${Barcode}/${Episode}.vcf | sort -k1,1 -k2,2n > bam_pass/${Barcode}/${Episode}_sorted.vcf
        grep "^#" bam_pass/${Barcode}/${Episode}.vcf > bam_pass/${Barcode}/${Episode}_header.vcf
        cat bam_pass/${Barcode}/${Episode}_header.vcf bam_pass/${Barcode}/${Episode}_sorted.vcf > bam_pass/${Barcode}/${Episode}.vcf
        rm bam_pass/${Barcode}/${Episode}_header.vcf bam_pass/${Barcode}/${Episode}_sorted.vcf

        # Merge BAM files
        samtools merge -@ 20 -u - bam_pass/${Barcode}/*bam | samtools sort -@ 20 -o bam_pass/${Barcode}/temp.bam && 
        samtools index bam_pass/${Barcode}/temp.bam && 
        samtools view -@ 20 -b bam_pass/${Barcode}/temp.bam ${Coordinate} > bam_pass/${Barcode}/${Episode}.bam && 
        samtools index bam_pass/${Barcode}/${Episode}.bam && 
        rm bam_pass/${Barcode}/temp*
        
        # Create coordinate BED file for this barcode/episode
        echo -e "${Coordinate}" | tr ':' '\t' | sed 's/-/\t/g' > bam_pass/${Barcode}/${Episode}_coordinate.bed

        # Run epi2me variant calling
        cd bam_pass/${Barcode}
        nextflow run epi2me-labs/wf-human-variation \
            --bam ${Episode}.bam \
            --ref /Users/javadjamshidi/Desktop/Refs/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
            --bed ${Episode}_coordinate.bed \
            --sample_name ${Episode} \
            --threads 24 \
            --bam_min_coverage 50 \
            --snp \
            --phased \
            --override_basecaller_cfg dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
            --annotation false \
            --out_dir variant_calling_output \
            --downsample_coverage false



        # Change to the original directory
        cd ../..
        
        # Change directory to the specific barcode folder before running the script
        cd bam_pass/${Barcode}
        python ../../phasing_automation_basecalling.py ${Episode}.bam ${Episode}.vcf
        # Return to the original directory
        cd ../..
    }
done < sample_sheet_nh