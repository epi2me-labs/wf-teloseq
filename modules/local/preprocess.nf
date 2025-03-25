/**
* Process reads, filtering out poor quality reads.workflow 
* Emits filtered reads, which span the telomere boundary.
**/
process process_reads {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    publishDir (
        "${params.out_dir}/${meta.alias}/unaligned_data/",
        mode: "copy",
    )

    input:
        tuple val(meta), path("reads.fastq"), path(stats)
    output:
        tuple val(meta), path("${processed_fastq}"), path("${summary_stats}"), path(stats)
    script:
        summary_stats = "${meta.alias}_telomere_unaligned_metrics.tsv";
        processed_fastq = "${meta.alias}_filtered_telomeric.fastq"
        """
        {
            samtools import -T '*' -OBAM -u reads.fastq || {
            echo '[ERROR] samtools import failed' >&2
            cat /dev/null
            }
        } \\
        | samtools reset -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl --no-PG -OBAM,level=1 \\
        | workflow-glue process_reads --summary-tsv-name ${summary_stats} ${meta.alias} - \\
        | samtools fastq -T '*' > ${processed_fastq}
        """
}
