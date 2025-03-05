/**
* Process reads, filtering out poor quality reads.workflow 
* Emits filtered reads, which span the telomere boundary.
**/
process process_reads {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    publishDir (
        "${params.out_dir}/${meta.alias}/",
        mode: "copy",
        overwrite: true
    )

    input:
        tuple val(meta), path("reads.fastq"), path(stats)
    output:
        tuple val(meta), path("reads_with_sub_telomere_no_error.fastq"), path("telomere_length_metrics.csv", optional: true), path("read_telomere_lengths.csv", optional: true), path(stats)
    script:
        """
        samtools import -T '*' -OBAM -u reads.fastq \\
        | samtools reset -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl --no-PG -OBAM,level=1 \\
        | workflow-glue process_reads \\
        | samtools fastq -T '*' > reads_with_sub_telomere_no_error.fastq
        """
}
