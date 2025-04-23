include { collapse_reference } from "./collapse_reference.nf"

params.reference_json_path = "input/reference.json"

workflow {
    def reference_json = file(params.reference_json_path)
    collapsed = collapse_reference(reference_json)
}
