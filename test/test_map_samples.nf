include { map_samples } from "../subworkflows/map_samples.nf"

params.fasta = "test/input/reps.fa"
params.samples = "test/input/samples/*"

workflow {
    map_samples(
        file(params.fasta),
        file(params.samples),
    )
}
