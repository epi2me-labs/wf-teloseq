/**
* Process reads, filtering out poor quality reads.workflow 
* Emits filtered reads, which span the telomere boundary.
**/
process process_reads {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq"), path(stats)
    output:
        tuple val(meta), path("reads_with_sub_telomere_no_error.fastq"), path("sample_raw_coverage.csv", optional: true), path("sample_raw_per_read_telomere_length.csv", optional: true), path(stats)
    script:
        """
        workflow-glue process_reads \\
            --input reads.fastq \\
            > reads_with_sub_telomere_no_error.fastq
        """
}
