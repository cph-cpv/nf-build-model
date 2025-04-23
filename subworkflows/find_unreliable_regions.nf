include { build_bowtie2_index } from "../modules/build_index.nf"
include { decompress } from "../modules/decompress.nf"

workflow find_unreliable_regions {
  take:
  host_fasta
  reference_fasta

  main:
  def find_mapped_regions_py = file("scripts/find_mapped_regions.py")
  def fragment_reference_py = file("scripts/fragment_reference.py")

  fragmented = fragment_reference(fragment_reference_py, reference_fasta)
  index = build_bowtie2_index(host_fasta)
  sam = run_fragment_mapping(fragmented, index)
  bam = sort_and_convert_sam(sam)
  unreliable_regions = find_mapped_regions(bam, find_mapped_regions_py)

  emit:
  unreliable_regions
}

process fragment_reference {
  conda 'env.yaml'
  cpus 1
  memory '500 MB'

  input:
  path fragment_reference_py
  path fasta

  output:
  path "fragmented.fa"

  script:
  """
  python3 ${fragment_reference_py} ${fasta} "fragmented.fa"
  """
}

process run_fragment_mapping {
  conda 'env.yaml'
  cpus 16
  memory '32 GB'

  input:
  path virus_segments
  path bowtie_index_path

  output:
  path "mapped.sam"

  script:
  """
  bowtie2 -x ${bowtie_index_path[0].baseName.replaceAll(/\.\d+/, '')} -p 15 -a -f ${virus_segments} -S mapped.sam
  """
}

process sort_and_convert_sam {
  debug true

  conda 'env.yaml'
  cpus 15
  memory '32 GB'

  input:
  path sam

  output:
  tuple path("mapped.bam"), path("mapped.bam.csi")

  script:
  """
  samtools sort --threads 15 --write-index -o mapped.bam ${sam}
  """
}

process find_mapped_regions {
  conda 'env.yaml'
  cpus 1
  publishDir "results/unreliable_regions", mode: 'copy'

  input:
  tuple path(bam), path(csi)
  path find_mapped_regions_py

  output:
  path "unreliable_regions.csv"

  script:
  """
  python3 ${find_mapped_regions_py} ${bam} "unreliable_regions.csv"
  """
}
