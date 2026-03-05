Barenbrug_NF_Pipeline

Automated Nextflow Pipeline for Variant Discovery and KASP Marker Development

This pipeline was designed by JASON CHARLES DAY as part of his MSc Thesis at The Centre for Bioinformatics and Computational Biology, Stellenbosch University, South Africa.

This pipeline is designed to discover high quality SNPs with the potential for use as KASP markers. The concept is split into two separate pipelines, both producing a set of sequences containing the relevant SNPs with surrounding sequence in a format that is suitable for KASP primer design submission as required by Barenbrug NZ. 

The process encapsulates and automates all steps of the variant calling and KASP marker identification process including:
-	Basic pre-processing of sequencing reads
-	Alignment of processed reads to a reference genome (Minimap2 for long reads and BWA MEM for short reads)
-	Sorting and Indexing of resulting BAM files using Samtools
-	Variant calling using Freebayes
-	Variant filtering using BCFtools
-	Extraction and generation of KASP markers + surrounding sequence using pysam
-	Standard end-run bioinformatics reporting using MultiQC
As previously stated, the pipeline is split into two primary streams, and this is dependant on the input type required by the user:
Stream 1: 
EVAK_pipe.nf Input: raw short- or long-read FASTQ 
Output: main/results/
Stream 2: 
EVAK_bam_pipe.nf Input: pre-aligned BAM files 
Output: main/results_from_bam/
DEPENDANCIES AND ENVIRONMENT
Core software
-	Nextflow >= 21.10.3
-	Docker default/latest version (both streams)
Software Tools (via Docker container images)
Stream 1 (EVAK_pipe.nf):
fastp:0.24.1
bwa:0.7.19
minimap2:2.29
samtools:1.21
mosdepth:0.3.10
freebayes:1.0.2
bcftools:1.21
pysam:0.22.0
multiqc:1.28
Stream 2 (EVAK_bam_pipe.nf):
samtools:1.21
mosdepth:0.3.10
freebayes:1.0.2
bcftools:1.21
pysam:0.22.0
multiqc:1.28
All are specified with explicit container images in EVAK_pipe.config and EVAK_bam_pipe.config.
PIPELINE STRUCTURE
STREAM 1 (EVAK_pipe.nf: raw FASTQ input)
1)	Read sample sheet.
2)	Pre-process short reads with fastp.
3)	Align reads: Short reads with BWA-MEM / Long reads with minimap2
4)	Sort and index BAM files with samtools.
5)	Calculate coverage using mosdepth.
6)	Call variants using FreeBayes.
7)	Filter SNPs using bcftools (SNP-only, quality metrics, depth, allele balance, allele frequency, flanking region).
8)	Retain only SNPs within annotated genes using a gene BED file.
9)	Extract KASP marker sequences using pysam.
10)	Aggregate QC reports with MultiQC.
STREAM 2 (EVAK_bam_pipe.nf: BAM input)
1)	Read BAM sample sheet.
2)	Index BAM files (samtools).
3)	Calculate coverage using mosdepth.
4)	Call variants using FreeBayes.
5)	Filter SNPs using bcftools (same filtering rules as Stream 1).
6)	Restrict SNPs to gene regions.
7)	Extract KASP marker sequences using pysam.
8)	Aggregate QC with MultiQC.
 
INPUTS
Reference genome FASTA (params.reference_genome)
Gene annotation BED file (params.genome_gene_bed_file)
Stream dependent sample sheet (samplesheet.csv / bam_samplesheet.csv)
These files to be placed in the main/data directory
STREAM 1 SAMPLE SHEET FORMAT (main/data/samplesheet.csv)
Required columns: id single_end read1 read2
Short-read example:
id,single_end,read1,read2 A2_test_NEA56_Contig1,false,/home/jasonday/Barenbrug_Nextflow_Pipeline/main/data/A2_test_R1.fastq.gz,/home/jasonday/Barenbrug_Nextflow_Pipeline/main/data/A2_test_R2.fastq.gz
Long-read example:
id,single_end,read1,read2 NEA6_long,true,/home/jasonday/Barenbrug_Nextflow_Pipeline/main/data/NEA6_long.fastq.gz,
(single_end=true indicates long reads; read2 remains empty)
STREAM 2 SAMPLE SHEET FORMAT (main/data/bam_samplesheet.csv)
Required columns: id bam_path
Example: 
id,bam_path NEA56_contig_7_bam_freebayes102,/home/jasonday/Barenbrug_Nextflow_Pipeline/main/data/ptg000007l_sorted.bam
 
OUTPUTS
STREAM 1 DEFAULT OUTPUT DIRECTORY: main/results/
Contains:
KASP marker outputs (main/results/kasp_markers)
MultiQC aggregated report (main/results/multiqc)
Nextflow pipeline_info (main/results/pipeline_info)
STREAM 2 DEFAULT OUTPUT DIRECTORY: main/results_from_bam/
Contains:
KASP marker outputs (main/results_from_bam/kasp_markers)
MultiQC report (main/results_from_bam/multiqc)
Nextflow pipeline info (main/results_from_bam/pipeline_info)
All intermediary results for both streams are written to main/work corresponding to unique Nextflow hash code.
 
RUNNING THE PIPELINES
From command line:
STREAM 1 (FASTQ INPUT)
Local run with Docker:
nextflow run EVAK_pipe.nf -c EVAK_pipe.config
STREAM 2 (BAM INPUT)
Local run with Docker:
nextflow run EVAK_bam_pipe.nf -c EVAK_bam_pipe.config
PARAMETER REFERENCE
--input Path to input samplesheet (FASTQ or BAM format)
--reference_genome Reference genome FASTA
--genome_gene_bed_file BED file with gene coordinates for gene-based filtering
--outdir Output directory
Variant filtering:
--min_allele_balance Minimum allele balance (default 0.2)
--max_allele_balance Maximum allele balance (default 0.8)
--min_allele_frequency Minimum allele frequency threshold (default 0.05)
--min_quality Minimum QUAL score (default 30)
--min_depth Minimum read depth (default 20)
--depth_multiplier Upper depth threshold = multiplier x mean depth
--flanking_region Number of flanking bases around SNP (default 100)
 
EXECUTION DEPENDANCIES
Default: local Docker execution
Execution profile frameworks exist for HPC-based runs and testing runs requiring minimal resources, however, need to be further fine-tuned prior to official use.
Docker is enabled by default in both streams and is required to be running on the system.
