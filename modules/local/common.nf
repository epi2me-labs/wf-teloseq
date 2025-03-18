// Force files to be collected before any consumers of the emitted channel try to use them
process collectFilesInDir {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("staging_dir/*"), val(dirname)
    output:
        tuple val(meta), path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}
