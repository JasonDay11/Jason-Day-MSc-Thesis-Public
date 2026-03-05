process FASTP {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_trimmed_{1,2}.fastq.gz"), emit: reads
    path "${meta.id}_fastp.json", emit: json
    path "${meta.id}_fastp.html", emit: html

    script:
    def prefix = "${meta.id}"
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${prefix}_trimmed_1.fastq.gz \
        -O ${prefix}_trimmed_2.fastq.gz \
        -j ${prefix}_fastp.json \
        -h ${prefix}_fastp.html
    """
}