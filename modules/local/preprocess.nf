/**
* Process reads, filtering out poor quality reads.workflow 
* Emits filtered reads, which span the telomere boundary.
**/
process process_reads {
    // NOTE This container has it's own EZcharts - remember to update!
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    publishDir (
        "${params.out_dir}/${meta.alias}/unaligned_data/",
        mode: "copy",
        pattern: "*.fastq"
    )
    publishDir (
        path: "${params.out_dir}/${meta.alias}/",
        mode: "copy",
        pattern: "*${stats_directory}*"
    )


    input:
        tuple val(meta), path("reads.fastq"), path("fastcat_stats")
    output:
        tuple val(meta), path("${processed_fastq}"), path("fastcat_stats"), path("${stats_directory}")
    script:
        summary_stats = "${meta.alias}_telomere_unaligned_metrics.tsv";  // nodef: used by output and glue
        processed_fastq = "${meta.alias}_filtered_telomeric.fastq";  // nodef: used by output and glue
        kde_stats = "${meta.alias}_kde_data.tsv";  // nodef: used by output and glue
        stats_directory = "stats";  // nodef: used to keep directory alignment between publishDir, glue and output
        """
        # default stats directory value
        mkdir stats
        {
            samtools import -T '*' -OBAM -u reads.fastq || {
            echo '[ERROR] samtools import failed' >&2
            cat /dev/null
            }
        } \\
        | samtools reset -x tp,cm,s1,s2,NM,MD,AS,SA,ms,nn,ts,cg,cs,dv,de,rl --no-PG -OBAM,level=1 \\
        | workflow-glue process_reads \\
            --summary-tsv-name ${summary_stats} \\
            --kde-tsv-name ${kde_stats} \\
            --base-stats-dir ${stats_directory} \\
            ${meta.alias} \\
            - \\
        | samtools fastq -T '*' > ${processed_fastq}
        """
}
