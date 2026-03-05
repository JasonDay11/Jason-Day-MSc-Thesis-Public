process SAMTOOLS_SORT_SHORT {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam
    path "${meta.id}.sorted.bam.stats", emit: stats

    script:
    """
    samtools sort -O bam -o ${meta.id}.sorted.bam ${sam}
    samtools stats ${meta.id}.sorted.bam > ${meta.id}.sorted.bam.stats
    """
}

process SAMTOOLS_INDEX_SHORT {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${bam}.bai"), emit: bai

    script:
    """
    samtools index ${bam}
    """
}