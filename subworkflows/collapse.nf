#!/usr/bin/env nextflow



workflow {
    def finish_py = file("finish.py")
    def organize_sequences_py = file("organize_sequences.py")
    def repair_py = file("repair.py")

    repaired_reference_path = repairReference(file("input/reference.json"), repair_py)

    otu_paths = organizeSequences(repaired_reference_path, organize_sequences_py) | flatten

    cluster_paths = otu_paths
        | flatMap { p ->
            p.listFiles().collect { fp -> tuple(p.baseName, fp) }
        }
        | clusterWithCdhit
        | collect

    finish(cluster_paths, repaired_reference_path, finish_py)
}

process repairReference {
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

process organizeSequences {
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

process clusterWithCdhit {
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
    debug true

    publishDir "results"

    input:
    path cluster_path
    path reference
    path finish_py

    output:
    path "reps.fa"
    path "reps_by_sequence.csv"
    path "summary.txt"

    script:
    """
    python3 ${finish_py} ${reference}
    """
}
