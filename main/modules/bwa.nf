process BWA_MEM {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.aligned.sam"), emit: sam

    script:
    def prefix = "${meta.id}"
    """
    # Index the reference genome if index doesn't exist
    if [ ! -f ${reference}.bwt ]; then
        bwa index ${reference}
    fi

    # Perform alignment
    bwa mem ${reference} ${reads[0]} ${reads[1]} > ${prefix}.aligned.sam
    """
}