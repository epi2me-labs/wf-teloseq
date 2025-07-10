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
        path: "${params.out_dir}/${meta.alias}/",
        mode: "copy",
        pattern: "*stats*"
    )

    input:
        tuple val(meta), path("reads.fastq"), path("reference.fasta"), path("stats")
    output:
        tuple val(meta), path("${output_bam}"), path("${output_bam}.csi"), emit: alignments
        tuple val(meta), path("stats"), emit: stats
    script:
        def summary_stats = "${meta.alias}_telomere_aligned_metrics.tsv";
        def boxplot_stats = "${meta.alias}_boxplot_values.tsv";
        def qc_stats = "${meta.alias}_qc_modes_metrics.tsv";
        def contig_stats = "${meta.alias}_contig_telomere_aligned_metrics.tsv";
        output_bam = "${meta.alias}_aligned_filtered_teloseqs.bam";  // nodef: used by output and glue
        """
        minimap2 -yax lr:hq --secondary=no -t ${task.cpus} --cap-kalloc 100m --cap-sw-mem 50m reference.fasta reads.fastq | samtools view -ubh | \\
            samtools sort -u - | \\
            workflow-glue process_alignments \\
                --summary-tsv-name ${summary_stats} \\
                --qc-tsv-name ${qc_stats} \\
                --boxplot-tsv-name ${boxplot_stats} \\
                --contig-summary-tsv-name ${contig_stats} \\
                --output-bam ${output_bam} \\
                --base-stats-dir stats \\
                ${meta.alias} \\
                -
        """
}
