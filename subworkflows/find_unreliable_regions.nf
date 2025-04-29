include { build_bowtie2_index } from "../modules/build_index.nf"
include { decompress } from "../modules/decompress.nf"

workflow find_unreliable_regions {
  take:
  host_fasta
  reference_fasta

  main:
  def combine_unreliable_regions_py = file("scripts/combine_unreliable_regions.py")
  def find_mapped_regions_py = file("scripts/find_mapped_regions.py")
  def find_unreliable_regions_nucleotide_percentage_r = file("scripts/find_unreliable_regions_nucleotide_percentage.r")
  def fragment_reference_py = file("scripts/fragment_reference.py")

  fragmented = fragment_reference(
    fragment_reference_py,
    reference_fasta,
  )

  check_fragments(
    fragmented,
    reference_fasta,
  )

  index = build_bowtie2_index(host_fasta)
  sam = run_fragment_mapping(fragmented, index)
  bam = sort_and_convert_sam(sam)

  mapped_unreliable_regions = find_mapped_regions(
    find_mapped_regions_py,
    bam,
  )

  nucleotide_percentage_unreliable_regions = find_unreliable_regions_nucleotide_percentage(
    find_unreliable_regions_nucleotide_percentage_r,
    reference_fasta,
  )

  unreliable_regions = combine_unreliable_regions(
    combine_unreliable_regions_py,
    mapped_unreliable_regions,
    nucleotide_percentage_unreliable_regions,
  )

  emit:
  unreliable_regions
}

process fragment_reference {
  conda 'env.yaml'
  cpus 1
  debug true
  memory '500 MB'
  publishDir "results/unreliable_regions"

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

process check_fragments {
  debug true

  input:
  path fragmented_fasta
  path reference_fasta

  script:
  """
  grep -c '^>' ${fragmented_fasta} | xargs echo "Fragment Count:"
  grep -c "^>" ${reference_fasta} | xargs echo "Input Fragment Count:"
  grep "^>" ${fragmented_fasta} | cut -f1 -d":" | uniq | wc -l | xargs echo "Fragmented Sequence Count:"
  """
}

process run_fragment_mapping {
  conda 'env.yaml'
  cpus 16
  debug true
  memory '32 GB'

  input:
  path fragments
  path bowtie_index_path

  output:
  path "mapped.sam"

  script:
  def index = bowtie_index_path[0].baseName.replaceAll(/\.\d+/, '')

  """
  bowtie2 -x ${index} -p 15 -a -f -U ${fragments} -S mapped.sam
  """
}

process sort_and_convert_sam {
  conda 'env.yaml'
  cpus 15
  memory '32 GB'
  publishDir "results/unreliable_regions"

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
  debug true
  publishDir "results/unreliable_regions"

  input:
  path find_mapped_regions_py
  tuple path(bam), path(csi)

  output:
  path "mapped_unreliable_regions.csv"

  script:
  """
  python3 ${find_mapped_regions_py} ${bam} "mapped_unreliable_regions.csv"
  """
}

process find_unreliable_regions_nucleotide_percentage {
  cpus 1
  memory '500 MB'

  input:
  path find_unreliable_regions_nucleotide_percentage_r
  path reference_fasta

  output:
  path "unreliable_regions_nucleotide_percentage.csv"

  script:
  """
  Rscript ${find_unreliable_regions_nucleotide_percentage_r} ${reference_fasta} unreliable_regions_nucleotide_percentage.csv
  """
}

process combine_unreliable_regions {
  cpus 1
  memory '500 MB'
  publishDir "results/unreliable_regions"

  input:
  path combine_unreliable_regions_py
  path unreliable_regions_mapped
  path unreliable_regions_nucleotide

  output:
  path "unreliable_regions.csv"

  script:
  """
  python3 ${combine_unreliable_regions_py} ${unreliable_regions_nucleotide} ${unreliable_regions_mapped} unreliable_regions.csv
  """
}
