process decompress {
    input:
    path y

    output:
    path { y.name.replace('.gz', '') }

    script:
    """
    export OUTPUT_NAME="${y.name.replace('.gz', '')}"

    if [[ ${y} == *.gz ]]; then
        gzip -df ${y} > \$OUTPUT_NAME
    fi
    """
}
