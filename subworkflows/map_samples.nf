include { build_bowtie2_index } from "../modules/build_index.nf"

workflow map_samples {
    take:
    fasta
    samples

    main:
    index = build_bowtie2_index(fasta)

    combined = Channel.from(samples) | combine(index) | run_sample_mapping 
    sorted = combined | sort_sample_mapping
}


process run_sample_mapping {
    cpus 15
    debug true
    publishDir "results/mapping"
    memory '40 GB'

    input:
    path files

    output:
    path "mapped.sam"

    script:
    """
    bowtie2 -x ${files[1].baseName.replaceAll(/\.\d+/, '')} -p 15 -a -U ${files[0]} -S mapped.sam
    """
}

process sort_sample_mapping {
    cpus 15
    memory '40 GB'

    input:
    path sam

    output:
    path "mapped.bam"

    script:
    """
    samtools view --threads 15 -b ${sam} > mapped.bam
    """
}
