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

  collapsed_fasta | view

  find_unreliable_regions(host, collapsed_fasta)
  map_samples(collapsed_fasta, samples)
}
