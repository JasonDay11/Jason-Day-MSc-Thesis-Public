process MINIMAP2 {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("${meta.id}.aligned.sam"), emit: sam

    script:
    """
    minimap2 -ax map-ont ${reference} ${reads} > ${meta.id}.aligned.sam
    """
}