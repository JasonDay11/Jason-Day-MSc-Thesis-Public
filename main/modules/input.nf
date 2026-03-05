process INPUT {
    input:
    path samplesheet

    output:
    path 'samples.csv'

    script:
    """
    awk -F ',' '{print \$1","\$2","\$3","\$4}' ${samplesheet} > samples.csv
    echo "Debug: Contents of samples.csv:"
    cat samples.csv
    """
}