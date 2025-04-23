include { decompress } from "../modules/decompress.nf"

workflow collapse_reference {
    take:
    reference_json_path

    main:
    def finish_py = file("scripts/finish.py")
    def organize_sequences_py = file("scripts/organize_sequences.py")
    def repair_py = file("scripts/repair.py")

    repaired_reference_path = repair_reference(
        reference_json_path,
        repair_py,
    )

    otu_paths = organize_sequences(repaired_reference_path, organize_sequences_py) | flatten

    cluster_paths = otu_paths
        | flatMap { p ->
            p.listFiles().collect { fp -> tuple(p.baseName, fp) }
        }
        | cluster_with_cdhit
        | collect

    collapsed = finish(cluster_paths, repaired_reference_path, finish_py)

    emit:
    collapsed
}

process repair_reference {
    cpus 1
    memory "5 GB"

    input:
    path reference
    path repair_py

    output:
    path "reference.json"

    script:
    """
    python3 ${repair_py} ${reference}
    mv reference_repaired.json reference.json
    """
}

process organize_sequences {
    cpus 1
    memory "200 MB"

    input:
    path reference
    path organize_sequences_py

    output:
    path "output/*"

    script:
    """
    python3 ${organize_sequences_py} ${reference} output
    """
}

process cluster_with_cdhit {
    cache "lenient"
    cpus 3
    memory "12 GB"

    input:
    tuple val(otu_id), path(segment_path)

    output:
    path "output_*"

    script:
    def output_path = "output_${otu_id}_${segment_path.baseName}"
    """                          
    mkdir -p '${output_path}'    
    cd-hit-est -T 2 -d 0 -c 0.99 -M 10000 -i '${segment_path / "sequences.fa"}' -o '${output_path}/clustered.fa'
    """
}


process finish {
    cpus 1
    memory "5 GB"
    publishDir "results/collapse_reference"

    input:
    path cluster_path
    path reference
    path finish_py

    output:
    tuple path("reps.fa"), path("reps_by_sequence.csv"), path("summary.txt")

    script:
    """
    python3 ${finish_py} ${reference}
    """
}
