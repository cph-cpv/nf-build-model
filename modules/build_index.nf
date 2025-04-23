process build_bowtie2_index {
    cpus 4
    memory '8 GB'

    input:
    path fasta_path

    output:
    path '*.bt2*'

    script:
    """
    bowtie2-build --threads 6 ${fasta_path} ${fasta_path.baseName}
    """
}
