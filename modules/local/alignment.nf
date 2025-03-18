process align_and_process {
    label "wf_teloseq"
    cpus params.alignment_threads
    memory "7 GB"
    publishDir (
        path: "${params.out_dir}/${meta.alias}/aligned_data/",
        mode: "copy",
        pattern: "*.bam*"
    )
    publishDir (
        path: "${params.out_dir}/${meta.alias}/aligned_data/",
        mode: "copy",
        pattern: "*_metrics.tsv"
    )

    input:
        tuple val(meta), path("reads.fastq"), path("reference.fasta"), path(stats)
    output:
        tuple val(meta), path("alignment_filtered_teloseqs.bam"), path("alignment_filtered_teloseqs.bam.csi"), emit: alignments
        tuple val(meta), path(stats), path("${summary_stats}"), path("${boxplot_stats}"), path("${qc_stats}"), path("${contig_stats}"), emit: alignment_stats
    script:
        summary_stats = "${meta.alias}_telomere_aligned_metrics.tsv";
        boxplot_stats = "${meta.alias}_boxplot_values.tsv";
        qc_stats = "${meta.alias}_qc_modes_metrics.tsv";
        contig_stats = "${meta.alias}_contig_telomere_aligned_metrics.tsv";
        """
        minimap2 -yax lr:hq --secondary=no -t ${task.cpus} --cap-kalloc 100m --cap-sw-mem 50m reference.fasta reads.fastq | samtools view -ubh | \\
            samtools sort -u - | \\
            workflow-glue process_alignments \\
                --summary-tsv-name ${summary_stats} \\
                --qc-tsv-name ${qc_stats} \\
                --boxplot-tsv-name ${boxplot_stats} \\
                --contig-summary-tsv-name ${contig_stats} \\
                --output-bam alignment_filtered_teloseqs.bam \\
                ${meta.alias} \\
                -
        """
}
