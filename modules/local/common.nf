process mapping {
    label "wf_teloseq"
    cpus 6
    memory "2.GB"

    input:
        tuple val(meta), path("reads.fastq"), path("reference.fasta")

    output:
        tuple val(meta), path("telomere.q${params.mapq}.bam"), path("telomere.q${params.mapq}.bam.bai"), emit: alignment
        tuple val(meta), path("reference.fasta"), emit: mappingref

    script:
        """
        minimap2 -ax map-ont -t ${task.cpus} reference.fasta reads.fastq | samtools view -bhq ${params.mapq} - | samtools sort -o telomere.q${params.mapq}.bam && samtools index telomere.q${params.mapq}.bam
        """
}
process fastq_stats {
    label "wf_teloseq"
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}/${meta.alias}/read_stats", mode: 'copy', overwrite: true
    input:
        tuple val(meta), path("reads.fastq")
        val output_name
        val header

    output:
        tuple val(meta), path(output_name)

    script:
        """
        seqkit stats -a reads.fastq \
        | awk 'NR==2 { \$1="${header}" }1' \
        | awk 'BEGIN {OFS=" "} { \$2=\$3=\$12=""; print }' \
        | tr -s ' ' '\t' > ${output_name}
        """
}

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
