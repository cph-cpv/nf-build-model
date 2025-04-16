include { buildBowtie2Index} from "../modules/build_index.nf"

workflow map_samples {
    take:
    collapsed_fasta
    samples

    main:
    collapsed_index = buildBowtie2Index(collapsed_fasta)
    mapping = run_mapping(samples, collapsed_index)

    row = Channel
        .fromPath("data/sample_name_paths_test.csv")
        .splitCsv()

    mapping = run_mapping(row, collapsed_index)

    emit:
    mapping.out
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
