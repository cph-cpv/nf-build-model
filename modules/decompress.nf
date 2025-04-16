process decompress {
    input:
    path y

    output:
    path y.baseName

    script:
    """
    if [[ ${y} == *.gz ]]; then
        gzip -d ${y}
    else
        cp ${y} ${y.baseName}
    fi
    """
}