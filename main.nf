include { build_bowtie2_index } from './modules/build_index.nf'
include { decompress as decompress_host } from "./modules/decompress.nf"
include { decompress as decompress_reference } from "./modules/decompress.nf"

include { find_unreliable_regions } from './subworkflows/find_unreliable_regions.nf'
include { collapse_reference } from './subworkflows/collapse_reference.nf'
include { map_samples } from './subworkflows/map_samples.nf'

params.host = "input/arabidopsis_thaliana.fa.gz"
params.reference = "input/reference.json"
params.samples = "input/samples/*"

workflow {
  def host = decompress_host(file(params.host))
  def reference = decompress_reference(file(params.reference))
  def samples = file(params.samples)

  collapsed = collapse_reference(reference)
  collapsed_fasta = collapsed | map { it[0] }

  extract_nucleotide_info(collapsed_fasta)
  find_unreliable_regions(host, collapsed_fasta)
  map_samples(collapsed_fasta, samples)
}

process extract_nucleotide_info {
  publishDir "results/nucleotide_info"

  input:
  path reference_json_path

  output:
  path "nucleotide_info.csv"

  script:
  """
  python3 scripts/extract_sequence_info.py ${reference_json_path}
  """

}

process run_mapping {
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

