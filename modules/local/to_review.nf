
process check_reference {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reference")  // can be compressed or uncompressed
    output:
        tuple val(meta), path("reference.fasta"), emit: reference
    script:
    """
    ##Extending telomere to each contig of the reference, as it has been observed misclassification of primary and secondary alignments
    ##if some arms have telomeres longer than the read but other contigs don't and are similar sequence but the correct site,
    ##it will map to the incorrect arm because it has a longer alignment from the prescence of telomere causing
    ##an expected primary alignment to become secondary and vice versa.

    workflow-glue extract_reference \\
        reference \\
        TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC \\
        ${params.restriction_site} \\
        ${params.beyond_cut} \\
        > extracted.fasta

    workflow-glue extend_telomere extracted.fasta TAACCC ${params.telomere_extension} > reference.fasta
    """
}



process trim_restriction_site {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("reads_restrict_trimmed.fastq"), emit: fastqrestrimmed
    script:
    """
    cutadapt -a "${params.restriction_site}" -e 0 -o reads_restrict_trimmed.fastq reads.fastq
    """
}

process trim_restriction_site_ref {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reference.fasta")
    output:
        tuple val(meta), path("reference_trim.fasta"), emit: reference
        tuple val(meta), path("cut_first.bed"), emit: cutsites
    script:
        """
        cutadapt -a "${params.restriction_site}" -e 0 -o reference_trim.fasta reference.fasta
        samtools faidx reference_trim.fasta
        awk -F'\t' '{print \$1"\t0\t"\$2}' reference_trim.fasta.fai > cut_first.bed
        """
}


process coverage_filter2 {
    label "wf_teloseq"
    cpus   4
    memory '2 GB'
    input:
        tuple val(meta), path("input.bam"), path("input.bam.bai"), path("input_reference.fasta"), val(coverage)
    output:
        tuple val(meta), path("./filtered/filtered.bam"),path("./filtered/filtered.bam.bai"), emit: alignment
        tuple val(meta), path("./filtered/filtered_reference.fasta"), emit: reference
    //remove 2nd and third column not needed. and rename first column for report
    script:
        """ 
        workflow-glue filter_bam_reference input.bam input_reference.fasta filtered $coverage $task.cpus
        """
}

process nm_filter {
    label "wf_teloseq"
    cpus 2
    memory '2 GB'
    input:
        tuple val(meta), path("telomere.bam"),path("telomere.bam.bai")
    output:
        tuple val(meta), path("telomere_raw.bam"),path("telomere_raw.bam.bai"), emit: alignment
    script:
    """
    workflow-glue filter_nm_reads telomere.bam telomere_raw.bam 4.000
    samtools index telomere_raw.bam
    """
}


//identify telomere end in reference
//get last telomere site in ref
//tidy up file and remove header
process telomere_sites{
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("tel_reference.fasta")
    output:
        tuple val(meta), path("telomere_boundaries.bed"), emit: telomerebed
    script:
        """
        workflow-glue get_telomere_boundaries_in_fastx tel_reference.fasta \
            > telomere_boundaries.bed
        """
}

process filtering {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai"), path("cut_sites.bed"), path("telomere_sites.txt")
        val(telomere_margin)
    output:
        tuple val(meta), path("${meta.alias}.tagged.bam"), path("${meta.alias}.tagged.bam.bai"), emit: finalbam

    script:
    //Identify reads to be removed from bam based upon filtering criteria.
    //['chr21_PATERNAL_P', 'chr21_MATERNAL_P'] have very long cut sitesm so 45000 length cut site limit is set so then strict filter does not just remove these type of contigs reads.
    //no filter is no reads removed 
    //low filter is keep only reads in which there end mapping position is 2000 bp beyond last telomere motif unless restriction site is before this
    //high filter is keep only reads in which there start mapping position is before last telomere motif identification and end is within 25 bp of cutsite
    """
    seqkit bam align.bam 2> sk-bam.tsv

    workflow-glue filter_bam_reads_output_id \\
        sk-bam.tsv \\
        cut_sites.bed \\
        telomere_sites.txt \\
        $telomere_margin \\
        25 \\
        high-keep-ids.txt \\
        low-keep-ids.txt \\
        no-extra-filter-keep-ids.txt

    # Convert filtering outputs to QC tags
    awk '{print \$1 "\\tHS"}' high-keep-ids.txt > high-tags.tsv
    awk '{print \$1 "\\tLS"}' low-keep-ids.txt > low-tags.tsv
    awk '{print \$1 "\\tNS"}' no-extra-filter-keep-ids.txt > no-filter-tags.tsv

    workflow-glue qc_tags align.bam tagged.bam high-tags.tsv low-tags.tsv no-filter-tags.tsv

    #Rename tagged BAM files to include the sample alias
    mv tagged.bam ${meta.alias}.tagged.bam
    mv tagged.bam.bai ${meta.alias}.tagged.bam.bai

    #Index the newly named tagged BAM file
    samtools index ${meta.alias}.tagged.bam
    """
}

process results {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta),
            path("tagged.bam"),
            path("tagged.bam.bai"),
            path("raw_coverage.csv"),
            val(coverage),
            path("mapping_ref.fasta")
    output:
        tuple val(meta),
            path("results_HS_chr_arm_coverage.csv"),
            path("results_HS_per_read_telomere_length.csv"),
            path("results_LS_chr_arm_coverage.csv"),
            path("results_LS_per_read_telomere_length.csv"),
            path("results_NS_chr_arm_coverage.csv"),
            path("results_NS_per_read_telomere_length.csv"),
            path("results_combined_summary.csv"),
            emit: for_report
    script:
    //search for locations of telomere sequences (x5 repeats) in individual reads.
    //Reverse locations file to select last occurance of each telomere match, thereby selecting end position of telomere.
    """
    samtools faidx mapping_ref.fasta

    seqkit bam tagged.bam -f Read,Ref,Acc,Strand,IsSec,IsSup \\
        2> tagged.sk-bam.tsv
    samtools fastq -T ZQ tagged.bam > tagged.fastq 
    workflow-glue get_telomere_boundaries_in_fastx tagged.fastq > telomere_read_length.txt
    workflow-glue summarise_telomere_lengths \\
        --mode alignment \\
        --unaligned raw_coverage.csv \\
        --seqkit_bam tagged.sk-bam.tsv \\
        --ref_fai mapping_ref.fasta.fai \\
        --min_coverage $coverage \\
        --telomere_lengths telomere_read_length.txt \\
        --output_prefix results
    """
}