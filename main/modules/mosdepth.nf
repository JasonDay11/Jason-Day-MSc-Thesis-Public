process MOSDEPTH {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.mosdepth.global.dist.txt"), emit: global
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), emit: summary
    path "${meta.id}.mosdepth.region.dist.txt"

    script:
    """
    mosdepth -n --fast-mode --by 500 ${meta.id} ${bam}
    """
}