#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SAMTOOLS_INDEX } from './modules/samtools'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { FREEBAYES } from './modules/freebayes.nf'
include { BCFTOOLS_FILTER_SNP; BCFTOOLS_FILTER_QUALITY } from './modules/bcftools.nf'
include { BCFTOOLS_FILTER_GENES } from './modules/gene_filter.nf'
include { EXTRACT_KASP_MARKERS } from './modules/kasp_markers.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow {
    // Create channel for input BAM files from samplesheet
    channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.id]
            def bam_file = file(row.bam_path)
            [meta, bam_file]
        }
        .set { ch_bam_input }

    // Index BAM files (in case they're not already indexed)
    SAMTOOLS_INDEX(ch_bam_input)
    
    // Join BAM and BAI files
    ch_bam_bai = ch_bam_input.join(SAMTOOLS_INDEX.out.bai, by: [0])

    // Run Mosdepth for coverage analysis
    MOSDEPTH(ch_bam_bai)

    // Call variants using Freebayes
    FREEBAYES(ch_bam_bai, params.reference_genome)

    // Filter variants for SNPs
    BCFTOOLS_FILTER_SNP(FREEBAYES.out.vcf)

    // Filter SNPs by quality metrics
    ch_vcf_mosdepth = BCFTOOLS_FILTER_SNP.out.vcf.join(MOSDEPTH.out.summary)
    BCFTOOLS_FILTER_QUALITY(
        ch_vcf_mosdepth,
        params.min_allele_balance,
        params.max_allele_balance,
        params.min_allele_frequency,
        params.min_quality,
        params.min_depth,
        params.depth_multiplier,
        params.flanking_region
    )

    // Filter SNPs to only include those within gene regions
    BCFTOOLS_FILTER_GENES(
        BCFTOOLS_FILTER_QUALITY.out.vcf,
        file(params.genome_gene_bed_file)
    )

    // Extract KASP markers from gene-filtered SNPs
    // Join gene-filtered VCF with original FREEBAYES VCF (contains ALL variants for masking)
    ch_kasp_input = BCFTOOLS_FILTER_GENES.out.vcf.join(FREEBAYES.out.vcf, by: [0])
    EXTRACT_KASP_MARKERS(
        ch_kasp_input,
        file(params.reference_genome),
        params.flanking_region
    )

    // Collect QC reports (excluding read processing steps)
    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary.flatten().filter(~/.*\.txt$/))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global.flatten().filter(~/.*\.txt$/))
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_FILTER_QUALITY.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_FILTER_GENES.out.stats)

    // Run MultiQC
    MULTIQC(ch_multiqc_files.collect())
}