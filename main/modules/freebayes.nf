process FREEBAYES {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf

    script:
    """
    freebayes -f ${reference} ${bam} > ${meta.id}.vcf
    """
}