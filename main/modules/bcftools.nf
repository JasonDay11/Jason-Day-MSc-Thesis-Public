process BCFTOOLS_FILTER_SNP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.snps.vcf"), emit: vcf

    script:
    """
    bcftools view -v snps ${vcf} > ${meta.id}.snps.vcf
    """
}

process BCFTOOLS_FILTER_QUALITY {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(vcf), path(mosdepth_summary)
    val min_allele_balance
    val max_allele_balance
    val min_allele_frequency
    val min_quality
    val min_depth
    val depth_multiplier
    val flanking_region

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf"), emit: vcf
    path "${meta.id}.filtered.vcf.stats", emit: stats

    script:
    """
    # Extract mean depth from the 'total' row
    mean_depth=\$(awk '\$1 == "total" {print \$4}' ${mosdepth_summary})
    if [ -z "\$mean_depth" ]; then
        echo "Error: Unable to extract mean depth from mosdepth summary file"
        exit 1
    fi
    echo "Extracted mean depth: \$mean_depth"
    
    # Calculate max depth using awk instead of bc
    max_depth=\$(awk -v mean=\$mean_depth -v mult=${depth_multiplier} 'BEGIN {print int(mean * mult)}')
    echo "Calculated max depth: \$max_depth"
    
    # Apply initial filters
    bcftools filter -i "AB >= ${min_allele_balance} && AB <= ${max_allele_balance} && AF >= ${min_allele_frequency} && QUAL >= ${min_quality} && FMT/DP >= ${min_depth} && FMT/DP <= \$max_depth" ${vcf} > temp.vcf

    # Extract positions of all SNPs
    bcftools query -f '%CHROM\\t%POS\\n' temp.vcf | sort -k2,2n > snp_positions.txt

    # Filter SNPs with at least ${flanking_region} bp flanking region free of other SNPs
    awk '{
        if (NR == 1) {
            prev = \$2;
            next;
        }
        if (\$2 - prev >= ${flanking_region}) {
            print \$1, \$2;
        }
        prev = \$2;
    }' snp_positions.txt > filtered_snp_positions.txt

    # Create a new VCF with only the filtered SNPs
    bcftools view -H temp.vcf | awk 'NR==FNR {pos[\$2]; next} (\$2 in pos)' filtered_snp_positions.txt - | \
        cat <(bcftools view -h temp.vcf) - > ${meta.id}.filtered.vcf
    
    # Generate stats for MultiQC
    bcftools stats ${meta.id}.filtered.vcf > ${meta.id}.filtered.vcf.stats

    # Cleanup
    rm snp_positions.txt filtered_snp_positions.txt temp.vcf
    """
}