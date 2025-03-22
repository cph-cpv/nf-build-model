// Default parameter input
reference_json = file("input/reference.json.gz")

process decompressReferenceJson {
    publishDir "results/decompress"
    tag "$y"

    input:
    path y

    output:
    path 'reference.json'

    script:
    """
    stat $y > stat_${y}
    gzip -df $y > reference.json
    """
}

process convertReferenceJsonToFasta {
    publishDir "results/fasta"

    input:
    path y

    output:
    path 'reference.fa'

    script:
    """
    jq -r '.otus | .[] | .isolates | .[] | .sequences | .[] | ">"+._id,.sequence' ${y} > reference.fa
    """
}

process buildBowtieIndex {
    publishDir "results/bowtie"

    cpus 6
    memory '4 GB'

    input:
    path y

    output:
    path 'reference.*.bt2'

    script:
    """
    bowtie2-build ${y} reference
    """
}

workflow {
  reference_json_path = decompressReferenceJson(reference_json)
  reference_fasta_path = convertReferenceJsonToFasta(reference_json_path)
  buildBowtieIndex(reference_fasta_path)
}