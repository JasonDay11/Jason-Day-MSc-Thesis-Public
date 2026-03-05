process BCFTOOLS_FILTER_GENES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(vcf)
    path gene_bed_file

    output:
    tuple val(meta), path("${meta.id}.gene_filtered.vcf"), emit: vcf
    path "${meta.id}.gene_filter.stats", emit: stats

    script:
    """
    # Count variants before filtering
    total_variants=\$(bcftools view -H ${vcf} | wc -l)
    echo "Total variants before gene filtering: \$total_variants"
    
    # Compress and index the VCF file first (required for -R option)
    bgzip -c ${vcf} > ${vcf}.gz
    bcftools index ${vcf}.gz
    
    # Filter VCF to only include variants within gene regions
    bcftools view -R ${gene_bed_file} ${vcf}.gz > ${meta.id}.gene_filtered.vcf
    
    # Count variants after filtering
    gene_variants=\$(bcftools view -H ${meta.id}.gene_filtered.vcf | wc -l)
    echo "Variants within gene regions: \$gene_variants"
    
    # Calculate percentage (fixed awk syntax)
    if [ \$total_variants -gt 0 ]; then
        percentage=\$(echo "\$gene_variants \$total_variants" | awk '{printf "%.2f", (\$1/\$2)*100}')
        echo "Percentage of variants in genes: \${percentage}%"
    else
        percentage="0.00"
        echo "No variants found for filtering"
    fi
    
    # Generate stats file
    cat > ${meta.id}.gene_filter.stats << EOF
# Gene Region Filtering Statistics for ${meta.id}
# Generated on: \$(date)
# BED file used: ${gene_bed_file}
Total_variants_before_filtering	\$total_variants
Variants_within_gene_regions	\$gene_variants
Percentage_in_genes	\${percentage}%
Variants_filtered_out	\$((\$total_variants - \$gene_variants))
EOF
    
    echo "Gene filtering completed for ${meta.id}"
    echo "Results: \$gene_variants/\$total_variants variants (\${percentage}%) are within gene regions"
    
    # Clean up temporary files
    rm -f ${vcf}.gz ${vcf}.gz.csi
    """
}