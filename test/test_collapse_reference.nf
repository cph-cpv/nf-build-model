include { collapse_reference } from "../subworkflows/collapse_reference.nf"
include { decompress as decompress_reference } from "../modules/decompress.nf"

params.reference_json = "test/input/reference.json.gz"

workflow {
    collapse_reference(
        decompress_reference(file(params.reference_json))
    )
}
