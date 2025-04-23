include { find_unreliable_regions } from "./find_unreliable_regions.nf"

params.host_fasta = "input/arabidopsis_thaliana.fa.gz"
params.reference_fasta = "input/reps.fa"

workflow {
    def host_fasta = file(params.host_fasta)
    def reference_fasta = file(params.reference_fasta)

    find_unreliable_regions(host_fasta, reference_fasta)
}
