segment_reference = file("scripts/segment_reference.py")
nucleotide_info = file("scripts/output/nucleotide_info.csv")
rscript_rle_conversion = file("scripts/rle_conversion.r")

process decompressHostFasta {
    publishDir "results/decompressHost"
    tag "$y"

    input:
    path y

    output:
    path 'host_genome.fa'

    script:
    """
    echo $y
    ls -l 
    stat $y > stat_${y}
    gzip -dfc ${y} > host_genome.fa
    """
}

process buildHostBowtieIndex {
    publishDir "results/hostIndex"

    cpus 6
    memory '4 GB'

    input:
    path host_fasta_path

    output:
    path 'host_genome_index.*.bt2'

    script:
    """
    bowtie2-build --threads 6 ${host_fasta_path} host_genome_index
    """
}

process runHostMapping {
    publishDir "results/hostMapping"

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
    
    ls -l

    samtools sort virus_segments_mapping.bam -o virus_segment_mapping.sorted.bam
    """ 
}

process featureExtractDataFrame{
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