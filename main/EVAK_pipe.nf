#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { INPUT } from './modules/input'
include { FASTP } from './modules/fastp'
include { BWA_MEM } from './modules/bwa'
include { MINIMAP2 } from './modules/minimap2'
include { SAMTOOLS_SORT; SAMTOOLS_INDEX } from './modules/samtools'
include { SAMTOOLS_SORT_SHORT; SAMTOOLS_INDEX_SHORT } from './modules/samtoolsshort'
include { SAMTOOLS_SORT_LONG; SAMTOOLS_INDEX_LONG } from './modules/samtoolslong'
include { MOSDEPTH } from './modules/mosdepth.nf'
include { FREEBAYES } from './modules/freebayes.nf'
include { BCFTOOLS_FILTER_SNP; BCFTOOLS_FILTER_QUALITY } from './modules/bcftools.nf'
include { BCFTOOLS_FILTER_GENES } from './modules/gene_filter.nf'
include { EXTRACT_KASP_MARKERS } from './modules/kasp_markers.nf'
include { MULTIQC } from './modules/multiqc.nf'

workflow {
    // Create channel for input samples
    INPUT(params.input)
    ch_samples = INPUT.out
        .splitCsv(header: true)
        .map { row -> 
            def meta = [id: row.id, single_end: row.single_end.toBoolean()]
            def reads = row.single_end.toBoolean() ? [file(row.read1)] : [file(row.read1), file(row.read2)]
            [meta, reads]
        }

    // Filter for short reads (paired-end) and long reads (single-end)
    ch_short_reads = ch_samples.filter { it -> !it[0].single_end }
    ch_long_reads = ch_samples.filter { it -> it[0].single_end }

    // Initialize empty fallback channels for BAM and BAI files
    // These will be overridden if reads are present
    ch_bam_short = channel.empty()
    ch_bai_short = channel.empty()
    ch_bam_long  = channel.empty()
    ch_bai_long  = channel.empty()


    // Process short reads: FASTP → BWA-MEM → sort → index
    // If no short reads are provided, the empty channel stays untouched
    FASTP(ch_short_reads)
    BWA_MEM(FASTP.out.reads, file(params.reference_genome))
    SAMTOOLS_SORT_SHORT(BWA_MEM.out.sam)
    SAMTOOLS_INDEX_SHORT(SAMTOOLS_SORT_SHORT.out.bam)

    ch_bam_short = SAMTOOLS_SORT_SHORT.out.bam
    ch_bai_short = SAMTOOLS_INDEX_SHORT.out.bai

    
    // Process long reads: Minimap2 → sort → index
    // If no long reads are provided, the empty channel stays untouched
    MINIMAP2(ch_long_reads, file(params.reference_genome))
    SAMTOOLS_SORT_LONG(MINIMAP2.out.sam)
    SAMTOOLS_INDEX_LONG(SAMTOOLS_SORT_LONG.out.bam)

    ch_bam_long = SAMTOOLS_SORT_LONG.out.bam
    ch_bai_long = SAMTOOLS_INDEX_LONG.out.bai


    // Combine both short and long BAM/BAI outputs using `mix()`
    // This ensures downstream steps receive unified channels even if one type is absent
    ch_bam = ch_bam_short.mix(ch_bam_long)
    ch_bai = ch_bai_short.mix(ch_bai_long)


    // Run Mosdepth
    ch_bam_bai = ch_bam.join(ch_bai, by: [0])
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

    // Collect all QC reports
    ch_multiqc_files = channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.flatten().filter(~/.*\.(json|html)$/))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary.flatten().filter(~/.*\.txt$/))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global.flatten().filter(~/.*\.txt$/))
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_SORT_SHORT.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_SORT_LONG.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_FILTER_QUALITY.out.stats)
    ch_multiqc_files = ch_multiqc_files.mix(BCFTOOLS_FILTER_GENES.out.stats)

    // Run MultiQC
    MULTIQC(ch_multiqc_files.collect())
}