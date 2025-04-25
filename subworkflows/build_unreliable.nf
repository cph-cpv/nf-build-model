include { buildBowtie2Index } from "../modules/build_index.nf"
include { decompress } from "../modules/decompress.nf"


workflow build_unreliable_regions {
  take:  
  host_fasta
  reference_fasta

  main:
  def fragment_reference_py = file("../scripts/fragment_reference.py")

  fragmented = fragment_reference(fragment_reference_py, reference_fasta)

  
  index = buildBowtie2Index(decompress(host_fasta))
  mapping = run_mapping(fragmented, index)

  row = Channel
    .fromPath('data/sample_name_paths_test.csv')
    .splitCsv()


  result = run_mapping(row, reference_index_path)

  


  rle_files = feature_extract_data_frame(
    result,
    reference_json_path,
    virus_segments,
    host_mapping
  )
}

process fragment_reference {
  input:
  path fragment_reference_py
  path fasta  

  output:
  path "fragmented.fa"

  script:
  """
  python ${fragment_reference_py} "fragmented.fa"
  """
}

process run_mapping {
    cpus 15
    memory '40 GB'

    input:
    path virus_segments
    path bowtie_index_path

    output:
    path "virus_segments_mapping.sorted.bam"
    path "virus_segments_mapping.bam"


    script:
    """
    bowtie2 -x ${bowtie_index_path[0].baseName.replaceAll(/\.\d+/, '')} -p 15 -a -f ${virus_segments} | samtools view -bS > virus_segments_mapping.bam
    samtools sort virus_segments_mapping.bam -o virus_segment_mapping.sorted.bam
    """ 
}




process feature_extract_data_frame{
  publishDir "results/rle_data"

  cpus 2

  input:
  path fileName
  path reference_json_path
  path virus_segments
  path host_mapping

  output:
  path "${fileName.baseName}_rle_data.rds"

  script:
  """
  pwd

  echo ${fileName.baseName}

  Rscript ${rscript_rle_conversion} ${fileName.baseName} ${reference_json_path} ${virus_segments} ${host_mapping} ${nucleotide_info}
  """
}

