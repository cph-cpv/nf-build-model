include { build_bowtie2_index } from "../modules/build_index.nf"

workflow map_samples {
    take:
    fasta
    samples

    main:
    index = build_bowtie2_index(fasta)

    Channel.from(samples) | combine(index) | run_sample_mapping
}


process run_sample_mapping {
    cpus 22
    debug true
    memory '40 GB'
    publishDir "results/mapping"

    input:
    path files

    output:
    path "*.bam*"

    script:
    def sample_name = files[0].baseName.replaceAll(/(\.fastq|\.fq)(\.gz)?$/, '')

    """
    bowtie2 -x ${files[1].baseName.replaceAll(/\.\d+/, '')} -p 16 -a -U ${files[0]} \
    | samtools view -bS - \
    | samtools sort -@ 4 -m 4G -o ${sample_name}.sorted.bam && samtools index ${sample_name}.sorted.bam
    echo ${sample_name}
    """
}
