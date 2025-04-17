include { buildBowtie2Index } from './modules/build_index.nf'
include { decompress } from "./modules/decompress.nf"

include { build_unreliable_regions } from './subworkflows/build_unreliable.nf' 
include { collapse_reference } from './subworkflows/collapse_reference.nf'
include { map_samples } from './subworkflows/map_samples.nf'

params.reference = "input/reference.json.gz"
params.samples_csv = "input/samples.csv"

workflow {
  def reference = decompress(file(params.reference))
  def samples = file(params.samples_path)

  collapsed = collapse_reference(reference)

  nucelotide_info = extract_nucleotide_info(collapsed)

  unreliable_regions = build_unreliable_regions(collapsed)

  sample_mapping = map_samples(row, collapsed_index)

  rle_samples = Channel.from(result)

  rle_files = featureExtractDataFrame(
    result,
    reference_json_path,
    virus_segments,
    host_mapping
  )
}

process extract_nucleotide_info {
  publishDir "results/nucleotideInfo"

  input:
  path reference_json_path

  output:
  path "nuceotide_info.csv"

  script:
  """
  python3 scripts/extract_nucleotide_info ${reference_json_path}
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

