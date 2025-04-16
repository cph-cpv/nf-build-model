// Default parameter input
include { decompressHostFasta; buildHostBowtieIndex; runHostMapping; featureExtractDataFrame } from './nextflow/build_unreliable.nf'
include { collapse_reference } from './subworkflows/collapse_reference.nf'


reference_json = file("input/reference.json.gz")
host_genome = file("input/arabidopsis_thaliana.fa.gz")
virus_segments = file("scripts/output/fragmented.fasta")


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

process buildReferenceBowtieIndex {
    publishDir "results/bowtie"

    cpus 6
    memory '4 GB'

    input:
    path y

    output:
    path 'reference.*.bt2'

    script:
    """
    bowtie2-build --threads 6 ${y} reference
    """
}

process runMapping {
    publishDir "results/mapping"

    cpus 15
    memory '40 GB'

    input:
    val row
    path bowtie_index_path

    output:
    path "${row[0]}.bam"


    script:
    """
    

    bowtie2 -x ${bowtie_index_path[0].baseName.replaceAll(/\.\d+/, '')} -p 15 -a -U ${row[1]} | samtools view -bS > ${row[0]}.bam
    """ 
}

workflow {
  // build unreliable_regions
  host_fasta = decompressHostFasta(host_genome)
  host_index_path = buildHostBowtieIndex(host_fasta)
  host_mapping = runHostMapping(virus_segments, host_index_path)


  reference_json_path = decompressReferenceJson(reference_json)
  reference_fasta_path = convertReferenceJsonToFasta(reference_json_path)
  bowtie_index_path = buildReferenceBowtieIndex(reference_fasta_path)
  

  row = Channel
    .fromPath('data/sample_name_paths_test.csv')
    .splitCsv()


  result  = runMapping(row, bowtie_index_path)

  rle_samples = Channel.from(result)

  rle_files = featureExtractDataFrame(result, reference_json_path, virus_segments, host_mapping)

}
