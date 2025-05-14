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
  def extract_sequence_info_py = file("scripts/extract_sequence_info.py")
  def match_sequence_info_py = file("scripts/match_sequence_info.py")

  def host = decompress_host(file(params.host))
  def labels = file(params.labels)
  def reference = decompress_reference(file(params.reference))
  def samples = file(params.samples)
  def samples_dir_path = samples[0].parent

  def associate_sample_labels_py = file("scripts/associate_sample_labels.py")
  sample_labels = associate_sample_labels(associate_sample_labels_py, labels, samples_dir_path)

  def extract_sample_viruses_py = file("scripts/extract_sample_viruses.py")
  sample_viruses = extract_sample_viruses(extract_sample_viruses_py, sample_labels, file(params.reference))

  collapsed = collapse_reference(reference)
  collapsed_fasta = collapsed | map { it[0] }

  sequence_info = extract_sequence_info(collapsed_fasta, extract_sequence_info_py)
  iimi_sequence_info = match_sequence_info_with_reference(match_sequence_info_py, sequence_info, file(params.reference))
  unreliable_regions = find_unreliable_regions(host, collapsed_fasta)
  mapped_samples = map_samples(collapsed_fasta, samples)

  rle_encode_mappings_r = file("scripts/rle_encode_mapping_data.r")
  rle_mappings = rle_encode_mappings(rle_encode_mappings_r, mapped_samples[0])

  rle_mappings_list = rle_mappings.collect()

  build_iimi_model_r = file("scripts/build_iimi_model.r")
  model = build_xgb_model(build_iimi_model_r, rle_mappings_list, unreliable_regions, iimi_sequence_info, sample_viruses)
}

process extract_sequence_info {
  conda 'env.yaml'
  publishDir "results/sequence_info"

  input:
  path reference_json_path
  path extract_sequence_info_py

  output:
  path "sequence_info.csv"

  script:
  """
  python3 ${extract_sequence_info_py} ${reference_json_path} sequence_info.csv
  """
}

process match_sequence_info_with_reference {
  conda 'env.yaml'
  publishDir "results/iimi_sequence_info"

  input:
  path match_sequence_info_py
  path sequence_info
  path reference_json_path

  output:
  path "iimi_sequence_info.csv"

  script:
  """
  python3 ${match_sequence_info_py} ${sequence_info} ${reference_json_path} iimi_sequence_info.csv
  """
}

/**
  * Associate user-provided sample labels with the sample names in the output files.
  */
process associate_sample_labels {
  cpus 1
  memory "5 GB"
  publishDir "results/sample_label_associations"

  input:
  path associate_sample_labels_py
  path labels
  path samples_dir

  output:
  path "sample_label_assocations.csv"

  script:
  """
  python3 ${associate_sample_labels_py} ${labels} ${samples_dir} "sample_label_assocations.csv"
  """
}

process extract_sample_viruses {
  cpus 1
  memory "5 GB"
  publishDir "results/sample_viruses"

  input:
  path extract_sample_viruses_py
  path sample_labels
  path reference

  output:
  path "sample_viruses.json"

  script:
  """
  python3 ${extract_sample_viruses_py} ${sample_labels} ${reference} "sample_viruses.json"
  """
}

process rle_encode_mappings {
  cpus 1
  memory "5 GB"
  publishDir "results/rle"

  input:
  path rle_encode_mappings_r
  path bam_file

  output:
  path "rle_${bam_file[0].getBaseName(2)}.rds"

  script:
  """
 Rscript ${rle_encode_mappings_r} ${bam_file[0]} .
 """
}

process build_xgb_model {
  publishDir "results/xgb_model"

  input:
  path build_iimi_model_r
  path rle_mappings
  path unreliable_regions
  path iimi_sequence_info
  path sample_viruses

  output:
  path "trained_xgb.rds"

  script:
  """
  sed -i '1 s/sequence_id/Virus segment/' ${unreliable_regions}
  Rscript ${build_iimi_model_r} . ${unreliable_regions} ${iimi_sequence_info} ${sample_viruses}
  """
}
