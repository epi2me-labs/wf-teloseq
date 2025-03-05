process results {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    publishDir "${params.out_dir}/${meta.alias}/read_stats", mode: 'copy', overwrite: true
    input:
        tuple val(meta),
            path("tagged.bam"),
            path("tagged.bam.csi"),
            path("sample_summary_metrics.csv"),
            path("mapping_ref.fasta")
    // Emit files named for each sample, as these are all dumped into the report
    output:
        tuple val(meta),
        path("${meta.alias}_results"),
        emit: for_report
    script:
    //search for locations of telomere sequences (x5 repeats) in individual reads.
    //Reverse locations file to select last occurance of each telomere match, thereby selecting end position of telomere.
    """
    seqkit bam tagged.bam -f Read,Ref,Acc,Strand,IsSec,IsSup \\
        2> tagged.sk-bam.tsv
    # Duplicated in extremis. I suppose the logic is these reads have now been trimmed?? 
    samtools fastq -T HT tagged.bam > tagged.fastq 
    workflow-glue get_telomere_boundaries_in_fastx tagged.fastq > telomere_read_length.txt
    workflow-glue summarise_telomere_lengths \\
        --unaligned sample_summary_metrics.csv \\
        --seqkit_bam tagged.sk-bam.tsv \\
        --min_coverage 5 \\
        --telomere_lengths telomere_read_length.txt \\
        --output_prefix ${meta.alias}

    mkdir -p ${meta.alias}_results
    mv *.csv ${meta.alias}_results
    """
}