include { find_unreliable_regions } from "../subworkflows/find_unreliable_regions.nf"

params.host_fasta = "test/input/arabidopsis_thaliana.fa.gz"
params.reference_fasta = "test/input/reps.fa"

workflow {
    find_unreliable_regions(
        file(params.host_fasta),
        file(params.reference_fasta)
    )
}
