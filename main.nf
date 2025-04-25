include { build_bowtie2_index } from './modules/build_index.nf'
include { decompress as decompress_host } from "./modules/decompress.nf"
include { decompress as decompress_reference } from "./modules/decompress.nf"

include { find_unreliable_regions } from './subworkflows/find_unreliable_regions.nf'
include { collapse_reference } from './subworkflows/collapse_reference.nf'
include { map_samples } from './subworkflows/map_samples.nf'

params.host = "input/arabidopsis_thaliana.fa.gz"
params.reference = "input/reference.json"
params.labels = "input/sample_labels.csv"
params.samples = "input/samples/*"

workflow {
  def host = decompress_host(file(params.host))
  def labels = file(params.labels)
  def reference = decompress_reference(file(params.reference))
  def samples = file(params.samples)
  def samples_dir_path = samples[0].parent

  def associate_sample_labels_py = file("scripts/associate_sample_labels.py")

  associate_sample_labels(associate_sample_labels_py, labels, samples_dir_path)

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

/**
  * Associate user-provided sample labels with the sample names in the output files.
  */
process associate_sample_labels {
  cpus 1
  debug true
  memory "5 GB"

  input:
  path associate_sample_labels_py
  path labels
  path samples_dir

  output:
  path "sample_label_assocations.csv"

  script:
  """
  ls
  python3 ${associate_sample_labels_py} ${labels} ${samples_dir} "sample_label_assocations.csv"
  """
}
