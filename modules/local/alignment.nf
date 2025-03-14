process align_and_process {
    label "wf_teloseq"
    cpus params.alignment_threads
    memory "7 GB"

    input:
        tuple val(meta), path("reads.fastq"), path("reference.fasta")

    output:
        tuple val(meta), path("alignment_filtered_teloseqs.bam"), path("alignment_filtered_teloseqs.bam.csi")

    script:
        """
        minimap2 -yax lr:hq -t ${task.cpus} --cap-kalloc 100m --cap-sw-mem 50m reference.fasta reads.fastq | samtools view -ubhq ${params.mapq} | \\
            samtools sort -u - | \\
            workflow-glue process_alignments --input-bam - \\
                --output-bam alignment_filtered_teloseqs.bam

        """
}
