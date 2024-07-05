#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include {
    getParams;
} from './lib/common'


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wf_teloseq"
    cpus 1
    memory 2.GB
    output: path "versions.txt"
    script:
    """
    python -c "import Bio; print(f'biopython,{Bio.__version__}')" >> versions.txt
    python -c "import matplotlib as mpl; print(f'matplotlib,{mpl.__version__}')" >> versions.txt
    python -c "import pyfastx; print(f'pyfastx,{pyfastx.__version__}')" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas as pd; print(f'pandas,{pd.__version__}')" >> versions.txt
    python -c "import numpy as np; print(f'numpy,{np.__version__}')" >> versions.txt
    seqkit version | grep "eqkit" | awk '{print "seqkit," \$2}' >> versions.txt
    python -c 'from importlib.metadata import version; v=version("edlib"); print(f"edlib,{v}")' \
        >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    csvtk version | awk '{print "csvtk," \$2}' >> versions.txt
    minimap2 --version | awk '{print "minimap2," \$1}' >> versions.txt
    cutadapt --version | awk '{print "cutadapt," \$1}' >> versions.txt
    vsearch --version 2>&1 | grep "vsearch " | sed 's/,.*//' | sed 's/ /,/' | sed 's/_.*//' >> versions.txt
    ( seqtk 2>&1 || true ) | grep "Version" | awk '{print "seqtk," \$2}' >> versions.txt
    """
}

//This process combines the first --denovo output reference with the user supplied contigs and error corrects and renames the entire combined reference.
process manualCuration {
    label "wf_teloseq"
    cpus = 6
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq"), path("coverage.txt")
        path("new_contigs.fasta")
        path("denovo_reference_round1.fasta")
        tuple val(meta2), path("reference_used_for_naming.fasta")
        val(enzyme)
    output:
        tuple val(meta), path("Manual_denovo_reference.fasta"), emit: ref1
    script:
    """

    extend_telomere.py new_contigs.fasta new_contigs2.fasta
    cutadapt --cut -6 denovo_reference_round1.fasta > reference2.fasta
    cat new_contigs2.fasta reference2.fasta > reference_manual.fasta
    minimap2 -ax map-ont -t $task.cpus reference_manual.fasta reads.fastq | samtools sort -o Telomere.bam
    samtools index Telomere.bam
    samtools faidx reference_manual.fasta
    freebayes -f reference_manual.fasta Telomere.bam >varnew.vcf
    bcftools view -i 'GT="1/1"' varnew.vcf -o filtered.vcf
    bgzip -c filtered.vcf > filtered.vcf.gz
    tabix -p vcf filtered.vcf.gz
    bcftools consensus -f reference_manual.fasta filtered.vcf.gz -o consensus.fasta

    minimap2 -ax map-ont -t $task.cpus consensus.fasta reads.fastq | samtools sort -o Telomere2.bam
    samtools index Telomere2.bam

    redundent_contig_removal_cov.py Telomere2.bam coverage.txt > idlisttoremove.txt
    #mosdepth -n -t 8 --fast-mode lowcov2 Telomere2.bam
    #awk -v cov="\$cov" '{if (\$4 > cov) print \$1}' lowcov2.mosdepth.summary.txt > idlisttokeep.txt
    #awk '{if (\$4 > cov) print \$1}' lowcov2.mosdepth.summary.txt > idlisttokeep.txt
    seqkit grep -v -f idlisttoremove.txt consensus.fasta > consensus_final.fasta

    #cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o correctedref.fasta reference_used_for_naming.fasta
    #blast_rename_ref.py -db correctedref.fasta -q consensus_final.fasta -o Manual_denovo_reference.fasta -i ${baseDir}/test_data/HG002_groupings.csv
    
    awk '/^>/{if(seq!=""){print seq "$enzyme"}; print; seq=""; next} {seq=seq""\$0} END{print seq "$enzyme"}' consensus_final.fasta > consensus_final2.fasta
    awk '/^>/ {print ">contig_"++i; next} {print}' consensus_final2.fasta > Manual_denovo_reference.fasta

    """
}

process coverage_calc {
    label "wf_teloseq"
    cpus 1
    memory 2.GB
    input:
        tuple val(meta), path(reads)
    output:
        tuple val(meta), path("coverage.txt")
    script:
    """
    if [[ "${params.mincoverage}" != -1 ]]; then
        echo "${params.mincoverage}" > coverage.txt
    else
        read_count=\$(wc -l < $reads | awk '{print \$1}')
        echo "\$((read_count / 4 / 92 * 20 / 100 ))" > coverage.txt
    fi
    """
}

process rmdup {
    label "wf_teloseq"
    cpus = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta), path("dedup.fastq")
    script:
    """
    seqkit rmdup -n reads.fastq.gz > dedup.fastq
    """
}

process rmdup2 {
    label "wf_teloseq"
    cpus = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq.gz")
    output:
        tuple val(meta), path("dedup.fastq")
    script:
    """
    seqkit rmdup -n reads.fastq.gz > dedup.fastq
    """
}

process subtelomere {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("subtelomere.txt")
    script:
    """
    cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o subtelomere.fastq reads.fastq
    awk 'NR%4 == 2 {print length(\$0)}' subtelomere.fastq > subtelomere2.txt
    tr -s ' ' '\t' < subtelomere2.txt > subtelomere.txt
    """
}

//remove short reads, could add quality filter too -Q, ALT pathway results in very short telomere so may not need this?
process remove_short1 {
    // TODO: ingress allows filtering of reads based on length + Q scores with fastcat
    // extra args
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("reads.short.fastq")
    script:
    """
    seqkit seq -m 100 -Q ${params.read_quality} reads.fastq -o reads.short.fastq
    """
}

//remove short reads, could add quality filter too -Q
process remove_short2 {
    // TODO: ingress allows filtering of reads based on length + Q scores with fastcat
    // extra args
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("reads.short.fastq")
    script:
    """
    seqkit seq -m 100 -Q ${params.read_quality} reads.fastq -o reads.short.fastq
    """
}

process check_reference {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        path(ref)
        val(enzyme)
    output:
        tuple val('reference_user'), path("reference.fasta"), emit: ref1
    publishDir "${params.out_dir}/reference/", mode: 'copy', overwrite: true
    // TODO: do we have to publish here instead of using the canonical way (with the
    // `output` process)?
    script:
    """
    ##Extending telomere to all references as observed misclassification of primary and secondary if some arms have telomeres
    ##long enough for mapping but the other similar site doesn't then that softclipping is taken into consideration
    ##and the primary can then become secondary and vice versa
    # TODO: the python scripts should probably use `pysam` instead of `biopython` as it
    # implicitly deals with compressed and uncompressed FASTx files
    if [[ $ref == *.gz ]]; then
        zcat $ref > reference.fasta 
        extract_reference.py reference.fasta $enzyme
        extend_telomere.py reference2.fasta reference.fasta
    elif [[ $ref == *.fasta ]]; then
        extract_reference.py $ref $enzyme
        extend_telomere.py reference2.fasta reference.fasta
    elif [[ $ref == *.fa ]]; then
        extract_reference.py $ref $enzyme
        extend_telomere.py reference2.fasta reference.fasta
    else
        pass
    fi
    """
}

process mapAndSeparateR1 {
    label "wf_teloseq"
    cpus = 6
    memory 2.GB
    input:
        tuple val(meta), path(input)
        tuple val(meta2), path("subtelomere_reference.fasta")
        val(cov)
    output:
        tuple val(meta), path('*.txt'), emit: clusters

    // to keep contigs to reference changed from this GGGGCGCGCAGCGCCGGCG
    script:
    """
    #remove telomere from reference used
    cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o Hg002refsubtelomere.fa subtelomere_reference.fasta
    #remove telomere telomere cleaned reads
    cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o trimmedreadscleaned.fastq $input
    #map cleaned subtelomere reads to subtelomere reference
    minimap2 -ax map-ont -t $task.cpus Hg002refsubtelomere.fa $input | samtools sort -o alignment.bam
    samtools index alignment.bam
    #get read information from mapped bam
    seqkit bam alignment.bam 2>alignmentseqkit
    #get fastq from mapped reads
    samtools fastq -F 4 alignment.bam > alignment.fastq
    #identify repeating motif
    seqkit locate --only-positive-strand -m 0 -p CGCCGGCGCAGGCG alignment.fastq > alignment.motif
    #summarise per read the repeating motif number
    tail -n +2 alignment.motif | cut -f1 | sort | uniq -c | awk '{print \$2 "\t" \$1}' > alignment.motifcounts
    #separate out mismapped reads, commented out left clip separation in this script as only useful in high coverage situations
    clustera.py alignmentseqkit alignment.motifcounts $cov
    """
}

process clusterAndExtractR1 {
    label "wf_teloseq"
    cpus = 4
    memory 2.GB
    input:
        tuple val(meta), path(clusterFiles), path(input)
        val(cov)
    output:
        tuple val(meta), path ("${clusterFiles}.set2.consensus_highestseqs.fasta") , emit: results

    script:
    """
    seqkit grep -f $clusterFiles $input > ${clusterFiles}.fastq
    seqtk seq -A ${clusterFiles}.fastq > ${clusterFiles}.set2.fasta
    vsearch --cluster_fast ${clusterFiles}.set2.fasta --strand plus --threads $task.cpus --maxseqlength 200000 --id 0.96 --consout ${clusterFiles}.set2.consensus.fasta
    extract_highest_seqs.py ${clusterFiles}.set2.consensus.fasta $cov

    # Check if output file is empty
    if [ ! -s ${clusterFiles}.set2.consensus_highestseqs.fasta ]; then
        # If empty, create an empty file
        touch ${clusterFiles}.set2.consensus_highestseqs.fasta
    fi
    """
}

process clusterAndExtractR2 {
    // TODO: instead of duplicating these processes, they could be defined in a separate
    // `.nf` file and imported several times
    label "wf_teloseq"
    cpus = 4
    memory 2.GB
    input:
        tuple val(meta), path(clusterFiles), path(input)
        val(cov)
    output:
        tuple val(meta), path ("${clusterFiles}.set2.consensus_highestseqs.fasta") , emit: results
    script:
    """
    seqkit grep -f $clusterFiles $input > ${clusterFiles}.fastq
    seqtk seq -A ${clusterFiles}.fastq > ${clusterFiles}.set2.fasta
    # --cluster_fast removed to reduce possibility of redundent contigs
    vsearch --cluster_fast ${clusterFiles}.set2.fasta --strand plus --threads $task.cpus --maxseqlength 200000 --id 0.96 --consout ${clusterFiles}.set2.consensus.fasta
    extract_highest_seqs.py ${clusterFiles}.set2.consensus.fasta $cov

    # Check if output file is empty
    if [ ! -s ${clusterFiles}.set2.consensus_highestseqs.fasta ]; then
        # If empty, create an empty file
        touch ${clusterFiles}.set2.consensus_highestseqs.fasta
    fi
    """
}

process clusterAndExtractR3 {
    label "wf_teloseq"
    cpus = 4
    memory 2.GB
    input:
        tuple val(meta), path(clusterFiles), path(input)
        val(cov)
        
    output:
        tuple val(meta), path ("${clusterFiles}.set2.consensus_highestseqs.fasta") , emit: results
    script:
    """
    seqkit grep -f $clusterFiles $input > ${clusterFiles}.fastq
    seqtk seq -A ${clusterFiles}.fastq > ${clusterFiles}.set2.fasta
    # --cluster_fast removed to reduce possibility of redundent contigs
    vsearch --cluster_fast ${clusterFiles}.set2.fasta --strand plus --threads $task.cpus --maxseqlength 200000 --id 0.96 --consout ${clusterFiles}.set2.consensus.fasta
    extract_highest_seqs.py ${clusterFiles}.set2.consensus.fasta $cov

    # Check if output file is empty
    if [ ! -s ${clusterFiles}.set2.consensus_highestseqs.fasta ]; then
        # If empty, create an empty file
        touch ${clusterFiles}.set2.consensus_highestseqs.fasta
    fi
    """
}

process combineRefR1 {
    label "wf_teloseq"
    cpus = 6
    memory 8.GB
    input: 
        tuple val(meta), path(input), path(resultFiles)
        val(cov)
    output:
        tuple val(meta), path("denovo_reference.fasta"), emit: ref1
        tuple val(meta), path("cutsites.bed"), emit: cutbed
        tuple val(meta), path("Telomere2.bam"), path("Telomere2.bam.bai"), emit: bam
    """
    cat $resultFiles | deduplicate.py new.fasta
    #trim all reads to TAACCC then add 5kbp of telomere
    extend_telomere.py new.fasta new2.fasta
    #map back original reads to polish 1|1 snp/indels
    minimap2 -ax map-ont -t $task.cpus new2.fasta $input | samtools sort -o Telomere.bam
    samtools index Telomere.bam
    samtools faidx new2.fasta
    freebayes -f new2.fasta Telomere.bam >varnew.vcf
    #script to filter but subset based upon coverage
    correct.py varnew.vcf new2.fasta predenovo_reference_temp.fasta
    cat predenovo_reference_temp.fasta | deduplicate.py predenovo_reference.fasta
    minimap2 -ax map-ont -t $task.cpus predenovo_reference.fasta $input | samtools sort -o Telomere2.bam
    samtools index Telomere2.bam

    redundent_contig_removal.py Telomere2.bam $cov > idlisttoremove.txt
    seqkit grep -v -f idlisttoremove.txt predenovo_reference.fasta > denovo_reference.fasta
    samtools faidx denovo_reference.fasta
    awk '{print \$1"\t"\$2"\t"\$2}' denovo_reference.fasta.fai > cutsites.bed
    """
}

process combineRefR2 {
    label "wf_teloseq"
    cpus = 6
    memory 8.GB
    input:
        tuple val(meta), path(input), path(resultFiles), path(ref), path("coverage.txt")
        tuple val(meta2), path("ref_original.fasta")
        val(enzyme)
    output:
        tuple val(meta), path("denovo_reference.fasta"), emit: ref1
        tuple val(meta), path("cutsites.bed"), emit: cutbed
        tuple val(meta), path("Raw_read.bam"), path("Raw_read.bam.bai"), emit: bam
    """
    cat $resultFiles > tmpref.fa
    extend_telomere.py tmpref.fa tmpref2.fa
    cat $ref tmpref2.fa > new2.fasta
    #map back original reads to polish 1|1 snp/indels
    minimap2 -ax map-ont -t $task.cpus new2.fasta $input | samtools sort -o Telomere.bam
    samtools index Telomere.bam
    samtools faidx new2.fasta
    #SHOULD I BE FILTERING BAM BEFORE THIS STAGE TO REMOVE LOW CONFIDENCE MAPPING?
    freebayes -f new2.fasta Telomere.bam >varnew.vcf
    #script to correct snp/indels in reference
    correct.py varnew.vcf new2.fasta predenovo_reference_temp.fasta
    #remove duplicate contigs using sequence that have come about by error in clustering contigs that already exist. Ignore first 80bp when comparing sequences
    cat predenovo_reference_temp.fasta | deduplicate.py predenovo_reference.fasta
    #map again
    minimap2 -ax map-ont -t $task.cpus predenovo_reference.fasta $input | samtools sort -o Raw_read.bam
    samtools index Raw_read.bam
    #remove duplicate contigs using coverage that have come about by error in clustering contigs that already exist
    redundent_contig_removal_cov.py Raw_read.bam coverage.txt > idlisttoremove.txt
    seqkit grep -v -f idlisttoremove.txt predenovo_reference.fasta > denovo_reference2.fasta

    #Add lost cut site for future pipeline use to end of contigs
    awk '/^>/{if(seq!=""){print seq "$enzyme"}; print; seq=""; next} {seq=seq""\$0} END{print seq "$enzyme"}' denovo_reference2.fasta > denovo_reference3.fasta
    #rename contigs ADD WITH CHR NAMING IN FUTURE
    awk '/^>/ {print ">contig_"++i; next} {print}' denovo_reference3.fasta > denovo_reference.fasta
    #create cutsites file showing where reads should map up to for later use
    samtools faidx denovo_reference.fasta
    awk '{print \$1"\t"\$2"\t"\$2}' denovo_reference.fasta.fai > cutsites.bed
    """
}

//combineRefR3
process combineRefR3 {
    label "wf_teloseq"
    cpus = 6
    memory 8.GB
    input: 
        tuple val(meta), path(input), path(resultFiles), path(ref), path("coverage.txt")
        tuple val(meta2), path("ref_original.fasta")
        val(enzyme)
    output:
        tuple val(meta), path("denovo_reference.fasta"), emit: ref1
        tuple val(meta), path("cutsites.bed"), emit: cutbed
        tuple val(meta), path("Raw_read.bam"), path("Raw_read.bam.bai"), emit: bam
    """
    cat $resultFiles > tmpref.fa
    extend_telomere.py tmpref.fa tmpref2.fa
    cat $ref tmpref2.fa > new2.fasta
    #map back original reads to polish 1|1 snp/indels
    minimap2 -ax map-ont -t $task.cpus new2.fasta $input | samtools sort -o Telomere.bam
    samtools index Telomere.bam
    samtools faidx new2.fasta

    freebayes -f new2.fasta Telomere.bam >varnew.vcf
    #script to correct snp/indels in reference
    correct.py varnew.vcf new2.fasta predenovo_reference_temp.fasta
    #remove duplicate contigs using sequence that have come about by error in clustering contigs that already exist. Ignore first 80bp when comparing sequences
    cat predenovo_reference_temp.fasta | deduplicate.py predenovo_reference.fasta
    #map again
    minimap2 -ax map-ont -t $task.cpus predenovo_reference.fasta $input | samtools sort -o Raw_read.bam
    samtools index Raw_read.bam
    #remove duplicate contigs using coverage that have come about by error in clustering contigs that already exist
    redundent_contig_removal_cov.py Raw_read.bam coverage.txt > idlisttoremove.txt
    seqkit grep -v -f idlisttoremove.txt predenovo_reference.fasta > denovo_reference2.fasta

    #rename final reference contigs, first remove variable regions #note could improve this further as can pair up within this naming and group further on similarity.
    #DONT USE AT MOMENT NEEDS FURTHER WORK
    #cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o correctedref.fasta ref_original.fasta
    #blast_rename_ref.py -db correctedref.fasta -q denovo_reference2.fasta -o denovo_reference5.fasta -i ${baseDir}/test_data/HG002_groupings.csv
    
    #Add lost cut site for future pipeline use to end of contigs
    awk '/^>/{if(seq!=""){print seq "$enzyme"}; print; seq=""; next} {seq=seq""\$0} END{print seq "$enzyme"}' denovo_reference2.fasta > denovo_reference3.fasta
    #rename contigs ADD WITH CHR NAMING IN FUTURE
    awk '/^>/ {print ">contig_"++i; next} {print}' denovo_reference3.fasta > denovo_reference.fasta
    #create cutsites file showing where reads should map up to for later use
    samtools faidx denovo_reference.fasta
    awk '{print \$1"\t"\$2"\t"\$2}' denovo_reference.fasta.fai > cutsites.bed
    """
}

process mapAndSeparateR2 {
    label "wf_teloseq"
    cpus = 6
    memory 2.GB
    input:
        tuple val(meta), path("alignment.bam"), path("alignment.bam.bai")
        val(cov)
    output:
        tuple val(meta), path('*_1.txt'), emit: clusters
    script:
    """
    #get read information from mapped bam
    seqkit bam alignment.bam 2>alignmentseqkit
    #get fastq from mapped reads
    samtools fastq alignment.bam > alignment.fastq
    #identify repeating motif but not using this now so could be removed with cluster script change...
    seqkit locate --only-positive-strand -m 0 -p CGCCGGCGCAGGCG alignment.fastq > alignment.motif
    #summarise per read the repeating motif number
    indel_count.py
    tail -n +2 alignment.motif | cut -f1 | sort | uniq -c | awk '{print \$2 "\t" \$1}' > alignment.motifcounts
    #separate out mismapped reads using just motif and pos
    cluster2b.py alignmentseqkit alignment.motifcounts $cov indel_counts.txt
    """
}

process mapAndSeparateR3 {
    label "wf_teloseq"
    cpus 2
    memory 2.GB
    input:
        tuple val(meta), path("alignment.bam"), path("alignment.bam.bai")
        val(cov)
    output:
        tuple val(meta), path('*_1.txt'), emit: clusters
    //set minimum 30 as clustered reads to produce reference and 20 coverage at the end to keep contigs to reference
    script:
    """
    #get read information from mapped bam
    seqkit bam alignment.bam 2>alignmentseqkit
    #get fastq from mapped reads
    samtools fastq alignment.bam > alignment.fastq
    #identify repeating motif but not using this now so could be removed with cluster script change...
    seqkit locate --only-positive-strand -m 0 -p CGCCGGCGCAGGCG alignment.fastq > alignment.motif
    #summarise per read the repeating motif number
    indel_count.py
    tail -n +2 alignment.motif | cut -f1 | sort | uniq -c | awk '{print \$2 "\t" \$1}' > alignment.motifcounts
    #separate out mismapped reads using just motif and pos
    cluster2b.py alignmentseqkit alignment.motifcounts $cov indel_counts.txt
    """
}

process trim_adapters {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("reads_trimmed.fastq"), emit: fastqtrimmed
    //trim adapter. why not use cutadapt, because adapter not basecalled well and precision important to end motifs
    """
    trim_adapter2.py reads.fastq reads_trimmed.fastq
    """

}

process filter_telomeres {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere.fastq")
    //identify telomere x10 repeat containing reads within the first 60-500bp as sequencing from 5' the telomere so should be there if telomere read
    """
    seqkit grep -s -R 60:500 -P -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" reads.fastq > telomere.fastq
    """
}

process filter_telomeres2 {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere.fastq")
    //identify telomere x10 repeat containing reads within the first 60-500bp as sequencing from 5' the telomere so should be there if telomere read
    """
    seqkit grep -s -R 60:500 -P -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" reads.fastq > telomere.fastq
    """
}

process reversecomplement {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere.fastq")
    //identify telomere x10 repeat containing reads within the first 60-500bp as sequencing from 5' the telomere so should be there if telomere read
    """
    seqtk seq -r reads.fastq > telomere.fastq
    """
}

process combinefastq {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
        tuple val(meta), path("reads2.fastq")
    output:
        tuple val(meta), path("telomere.fastq")
    //identify telomere x10 repeat containing reads within the first 60-500bp as sequencing from 5' the telomere so should be there if telomere read
    """
    cat reads.fastq reads2.fastq > telomere.fastq
    """
}

process filter_nontelomeres {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere_and_non_telomere.fastq")
    //filter further subset that should not have telomere sequence in the last 70bp as shortest cut site is beyond this so should have atleast this amount of sequence at end that
    //is not telomere and don't want telomere only reads as could be fragments or artefacts. 
    """
    seqkit grep -v -s -R -70:-1 -P -p "TAACCCTAACCCTAACCCTAACCC" reads.fastq > telomere_and_non_telomere.fastq
    """
}

process filter_motifs {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("removereadids.txt")
    //seqkit forward strand the telomere repeat and identify all locations for eachr read
    //reverse the file so last telomere location for each read is first and print one row for each read, effectively getting last telomere match position for each read
    //retrieve just id and length information I need
    //search for basecalling error in telomere
    //Idnetify reads with high intensity so lots of kmer error clustered error within the telomere as this will not map and get softclipped, plus not useful for telomere length.
    //ensure just one read ID and @ is removed as not needed later on.
    """
    seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC reads.fastq > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1" "\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt

    seqkit locate --only-positive-strand -m 0 -p GTATAG,CGCGCGCG,CCACCG,AGCGACAG,ATAAGT,CCTCGTCC,TATAGT,AGTACT,GAGTCC,TATAGT,TATACA,TGGTCC,CTCTCCTCT reads.fastq > motifs2.txt
    error_reads.py motifs2.txt locationstelomerelast2.txt errors.txt

    awk '!a[\$1]++' errors.txt > removereadids.txt
    sed -i 's/@//g' removereadids.txt

    """
}

//remove basecalling error reads from telomere only identified reads
process filter_motifs_reads1 {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq"), path("remove_ids.txt")
    output:
        tuple val(meta), path("Telomere_reads.fastq")
    """
    seqkit grep -v -f remove_ids.txt reads.fastq > Telomere_reads.fastq
    """
}

//remove basecalling error reads from telomere and subtelomere identified reads
process filter_motifs_reads2 {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq"), path("remove_ids.txt")
    output:
        tuple val(meta), path("Telomere_reads.fastq")
    """
    seqkit grep -v -f remove_ids.txt reads.fastq > Telomere_reads.fastq
    """
}

process fastq_stats {
    // TODO: could we use the stats produced by `fastcat` during ingress instead?
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("raw.txt")
    //remove 2nd and third column not needed. and rename first column for report
    """
    seqkit stats -a reads.fastq > stats3.txt
    awk 'NR==2 { \$1="Raw_Reads" }1' stats3.txt > stats2.txt
    awk 'BEGIN {OFS=" "} { \$2=\$3=\$12=""; print }' stats2.txt > raw2.txt
    tr -s ' ' '\t' < raw2.txt > raw.txt
    """
}

process fastq_stats2 {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere.txt")
    //remove 2nd and third column not needed. and rename first column for report
    """
    seqkit stats -a reads.fastq > stats3.txt
    awk 'NR==2 { \$1="Telomere_Reads" }1' stats3.txt > stats2.txt
    awk 'BEGIN {OFS=" "} { \$2=\$3=\$12=""; print }' stats2.txt > telomere2.txt
    tr -s ' ' '\t' < telomere2.txt > telomere.txt
    """
}

process fastq_stats3 {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere_subtelomere.txt")
    //remove 2nd and third column not needed. and rename first column for report
    """
    seqkit stats -a reads.fastq > stats3.txt
    awk 'NR==2 { \$1="Telomere_Subtelomere_Reads" }1' stats3.txt > stats2.txt
    awk 'BEGIN {OFS=" "} { \$2=\$3=\$12=""; print }' stats2.txt > telomere_subtelomere2.txt
    tr -s ' ' '\t' < telomere_subtelomere2.txt > telomere_subtelomere.txt
    """
}

process mappingbam {
    label "wf_teloseq"
    cpus = 6
    memory 2.GB
    input:
        tuple val(meta), path("reads.fastq")
        tuple val(meta2), path("mapping_reference2.fasta")
    output:
        tuple val(meta), path("telomere.q${params.mapq}.bam"),path("telomere.q${params.mapq}.bam.bai"), emit: alignments
        tuple val(meta), path("telomere.bam"), path("telomere.bam.bai"), emit: alignments2  // TODO: looks like this isn't used anywhere
        tuple val(meta), path("mapping_reference.fasta"), emit: mappingref
    script:
    """
    cp mapping_reference2.fasta mapping_reference.fasta  # TODO: why cp this?
    minimap2 -ax map-ont -t $task.cpus mapping_reference.fasta reads.fastq | samtools sort -o telomere.bam
    samtools index telomere.bam
    samtools view -bq ${params.mapq} -h telomere.bam > "telomere.q${params.mapq}.bam"
    samtools index "telomere.q${params.mapq}.bam"
    """
}

//identify enzyme cut sites in reference, change -p sequence if different enzyme used
//get first cut site in ref
//tidy up file and remove header
process cut_sites{
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta2), path("cut_reference.fasta")
        val(enzyme)
    output:
        tuple val(meta2), path("cutfirst.bed"), emit: cutbed
    script:
    """       
    seqkit locate --only-positive-strand -m 0 -p $enzyme cut_reference.fasta > cutsites2.txt
    awk '!a[\$1]++' cutsites2.txt > cutsites3.txt
    awk -F'\t' '{print \$1"\t"\$5"\t"\$5}' cutsites3.txt > cutfirst4.txt
    grep -v 'seqID' cutfirst4.txt > cutfirst.bed
    """
}

//identify telomere end in reference, change -p sequence if different enzyme used
//get lsat telomere site in ref
//tidy up file and remove header
process telomere_sites{
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta2), path("tel_reference.fasta")
    output:
        tuple val(meta2), path("locationstelomerelast2.bed"), emit: telomerebed
    script:
    """       
    seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC,AACCCTAACCCTAACCCTAACCCTAACCCT,ACCCTAACCCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCTAACCCTAA,CCTAACCCTAACCCTAACCCTAACCCTAAC,CTAACCCTAACCCTAACCCTAACCCTAACC tel_reference.fasta > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1"\t"\$5"\t"\$5}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > locationstelomerelast2.txt
    grep -v 'seqID' locationstelomerelast2.txt > locationstelomerelast2.bed
    """
}

process filtering {
    tag { meta.alias }
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
        tuple val(meta2), path("cutsites.txt")
        tuple val(meta2), path("telomeresites.txt")
    output:
        tuple val(meta), path("highfiltered.bam"), path("highfiltered.bam.bai"), emit: finalbam
        tuple val(meta), path("lowfiltered.bam"), path("lowfiltered.bam.bai"), emit: lowfinalbam
        tuple val(meta), path("nofiltered.bam"), path("nofiltered.bam.bai"), emit: nofinalbam
        tuple val(meta),
            path("highfiltered.bam"),
            path("highfiltered.bam.bai"),
            path("lowfiltered.bam"),
            path("lowfiltered.bam.bai"),
            path("nofiltered.bam"),
            path("nofiltered.bam.bai"),
            emit: combined
    script:
    //Identify reads to be removed from bam based upon filtering criteria.
    //['chr21_PATERNAL_P', 'chr21_MATERNAL_P'] have very long cut sitesm so 45000 length cut site limit is set so then strict filter does not just remove these type of contigs reads.
    //no filter is no reads removed 
    //low filter is keep only reads in which there end mapping position is 80 bp beyond last telomere motif
    //high filter is keep only reads in which there start mapping position is before last telomere motif identification and end is within 25 bp of cutsite
    """
    seqkit bam align.bam 2>seqkitbamout.txt
    filter_bam_reads_output_id.py --strict seqkitbamout.txt cutsites.txt telomeresites.txt id.txt

    awk '!seen[\$0]++' id.txt > id2.txt 
    samtools view -N id2.txt -o temp.bam align.bam 
    samtools index temp.bam
    samtools view -h -o highfiltered.bam temp.bam
    samtools index highfiltered.bam

    filter_bam_reads_output_id.py --lenient seqkitbamout.txt cutsites.txt telomeresites.txt lowid.txt 

    awk '!seen[\$0]++' lowid.txt > lowid2.txt 
    samtools view -N lowid2.txt -o lowtemp.bam align.bam
    samtools index lowtemp.bam
    samtools view -h -o lowfiltered.bam lowtemp.bam
    samtools index lowfiltered.bam

    filter_bam_reads_output_id.py --none seqkitbamout.txt cutsites.txt telomeresites.txt noneid3.txt

    awk '!seen[\$0]++' noneid3.txt > noneid4.txt 
    samtools view -N noneid4.txt -o nonetemp4.bam align.bam
    samtools index nonetemp4.bam
    samtools view -h -o nofiltered.bam nonetemp4.bam
    samtools index nofiltered.bam
    """
}

process raw_telomere_analysis {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("telomere.fastq")
    output:
        tuple val(meta), path("Sample_raw_Coverage.csv"), emit: covraw
        tuple val(meta), path("Sample_raw_Boxplot_of_Telomere_length.pdf"), emit: pdfraw
        tuple val(meta), path("Sample_raw_Per_Read_telomere_length.csv"), emit: plotraw
    script:
    //search for locations of telomere sequences (x5 repeats) in individual reads.
    //Reverse locations file to select last occurance of each telomere match, thereby selecting end position of telomere
    //subtract the first telomere location from last so remove adapter length from telomere length.
    """
    cat telomere.fastq | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC,AACCCTAACCCTAACCCTAACCCTAACCCT,ACCCTAACCCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCTAACCCTAA,CCTAACCCTAACCCTAACCCTAACCCTAAC,CTAACCCTAACCCTAACCCTAACCCTAACC > locationstelomere.txt
    tac locationstelomere.txt | awk '!a[\$1]++' > locationstelomerelast.txt
    awk -F'\t' '{print \$1" "\$6}' locationstelomerelast.txt | sort -r | tr ' ' '\t' > end
    awk -F'\t' '!a[\$1]++' locationstelomere.txt | awk -F'\t' '{print \$1" "\$5}' | sort -r | tr ' ' '\t' > start
    awk 'BEGIN {OFS="\t"} FNR==NR {if (NR>1) {a[\$1]=\$2}; next} FNR==1 {print} FNR>1 {\$2=\$2-a[\$1]; print}' start end > Sample
    telomerewindowV1.py telomere.fastq telomere_read_length.txt
    telomere_length_coverage_raw2.py Sample telomere_read_length.txt
    """
}

process results {
    tag { meta.alias }
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta),
            path("highfiltered.bam"),
            path("highfiltered.bam.bai"),
            path("lowfiltered.bam"),
            path("lowfiltered.bam.bai"),
            path("nofiltered.bam"),
            path("nofiltered.bam.bai"),
            path("raw_coverage.csv"),
            path("coverage.txt")
        tuple val(meta2), path("mapping_ref.fasta")
    output:
        tuple val(meta),
            path("highfiltered_chr_arm_Coverage.csv"),
            path("highfiltered_Per_Read_telomere_length.csv"),
            path("lowfiltered_chr_arm_Coverage.csv"),
            path("lowfiltered.csv"),
            path("lowfiltered_Per_Read_telomere_length.csv"),
            path("nofiltered_chr_arm_Coverage.csv"),
            path("nofiltered_Per_Read_telomere_length.csv"),
            path("output.csv"), 
            emit: for_report
        tuple val(meta),
            path("highfiltered_Boxplot_of_Telomere_length.pdf"),
            path("highfiltered_chr_arm_Boxplot_of_Telomere_length.pdf"),
            path("nofiltered_chr_arm_Boxplot_of_Telomere_length.pdf"),
            path("lowfiltered_chr_arm_Boxplot_of_Telomere_length.pdf"),
            path("nofiltered_Boxplot_of_Telomere_length.pdf"),
            path("lowfiltered_Boxplot_of_Telomere_length.pdf"),
            emit: pdf
        tuple val(meta),
            path("nofiltered_Per_Read_telomere_length.csv"),
            path("highfiltered_chr_arm_Coverage.csv"),
            path("lowfiltered_Per_Read_telomere_length.csv"),
            path("highfiltered_Per_Read_telomere_length.csv"),
            path("nofiltered_chr_arm_Coverage.csv"),
            path("lowfiltered_chr_arm_Coverage.csv"),
            emit: alldata
    script:
    //search for locations of telomere sequences (x5 repeats) in individual reads.
    //Reverse locations file to select last occurance of each telomere match, thereby selecting end position of telomere
    //remove the first location of telomere from the last to remove any adapter length contributing.
    """
    samtools index highfiltered.bam
    samtools fastq highfiltered.bam | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC,AACCCTAACCCTAACCCTAACCCTAACCCT,ACCCTAACCCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCTAACCCTAA,CCTAACCCTAACCCTAACCCTAACCCTAAC,CTAACCCTAACCCTAACCCTAACCCTAACC > locationstelomerestrict.txt
    tac locationstelomerestrict.txt | awk '!a[\$1]++' > locationstelomerelaststrict.txt
    awk -F'\t' '{print \$1" "\$6}' locationstelomerelaststrict.txt | sort -r | tr ' ' '\t' > strict_end
    awk -F'\t' '!a[\$1]++' locationstelomerestrict.txt | awk -F'\t' '{print \$1" "\$5}' | sort -r | tr ' ' '\t' > strict_start
    awk 'BEGIN {OFS="\t"} FNR==NR {if (NR>1) {a[\$1]=\$2}; next} FNR==1 {print} FNR>1 {\$2=\$2-a[\$1]; print}' strict_start strict_end > strict_tel_length
    samtools faidx mapping_ref.fasta
    seqkit bam highfiltered.bam 2>highfiltered
    samtools fastq highfiltered.bam > highfiltered.fastq
    telomerewindowV1.py highfiltered.fastq high_telomere_read_length.txt
    telomere_length_coverage3.py highfiltered strict_tel_length mapping_ref.fasta.fai coverage.txt high_telomere_read_length.txt

    samtools index lowfiltered.bam
    samtools fastq lowfiltered.bam | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC,AACCCTAACCCTAACCCTAACCCTAACCCT,ACCCTAACCCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCTAACCCTAA,CCTAACCCTAACCCTAACCCTAACCCTAAC,CTAACCCTAACCCTAACCCTAACCCTAACC > locationstelomerelenient.txt
    tac locationstelomerelenient.txt | awk '!a[\$1]++' > locationstelomerelastlenient.txt
    awk -F'\t' '{print \$1" "\$6}' locationstelomerelastlenient.txt | sort -r | tr ' ' '\t' > lenient_end
    awk -F'\t' '!a[\$1]++' locationstelomerelenient.txt | awk -F'\t' '{print \$1" "\$5}' | sort -r | tr ' ' '\t' > lenient_start
    awk 'BEGIN {OFS="\t"} FNR==NR {if (NR>1) {a[\$1]=\$2}; next} FNR==1 {print} FNR>1 {\$2=\$2-a[\$1]; print}' lenient_start lenient_end > lenient_tel_length
    seqkit bam lowfiltered.bam 2>lowfiltered
    samtools fastq lowfiltered.bam > lowfiltered.fastq
    telomerewindowV1.py lowfiltered.fastq low_telomere_read_length.txt 

    telomere_length_coverage3.py lowfiltered lenient_tel_length mapping_ref.fasta.fai coverage.txt low_telomere_read_length.txt

    samtools index nofiltered.bam
    samtools fastq nofiltered.bam | seqkit locate --only-positive-strand -m 1 -p TAACCCTAACCCTAACCCTAACCCTAACCC,AACCCTAACCCTAACCCTAACCCTAACCCT,ACCCTAACCCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCTAACCCTAA,CCTAACCCTAACCCTAACCCTAACCCTAAC,CTAACCCTAACCCTAACCCTAACCCTAACC > locationstelomereraw.txt
    tac locationstelomereraw.txt | awk '!a[\$1]++' > locationstelomerelastraw.txt
    awk -F'\t' '{print \$1" "\$6}' locationstelomerelastraw.txt | sort -r | tr ' ' '\t' > raw_end
    awk -F'\t' '!a[\$1]++' locationstelomereraw.txt | awk -F'\t' '{print \$1" "\$5}' | sort -r | tr ' ' '\t' > raw_start
    awk 'BEGIN {OFS="\t"} FNR==NR {if (NR>1) {a[\$1]=\$2}; next} FNR==1 {print} FNR>1 {\$2=\$2-a[\$1]; print}' raw_start raw_end > raw_tel_length
    seqkit bam nofiltered.bam 2>nofiltered
    samtools fastq nofiltered.bam > nofiltered.fastq
    telomerewindowV1.py nofiltered.fastq no_telomere_read_length.txt

    telomere_length_coverage3.py nofiltered raw_tel_length mapping_ref.fasta.fai coverage.txt no_telomere_read_length.txt

    combine_files.py 
    """
}

process makeReport {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        val metadata
        path "versions/*"
        path "params.json"
        path "data/*"
        val mappingreport
    output:
        path "wf-teloseq-*.html"
    script:
        String extra_arg = ""
        if (mappingreport) {
            extra_arg = "--mappingreport"
        }
        String report_name = "wf-teloseq-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json

    workflow-glue report $report_name \
        --metadata metadata.json \
        --versions versions \
        --params params.json \
        --data data \
        $extra_arg
    """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wf_teloseq"
    cpus 1
    memory 2.GB
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

process collectFilesInDir {
    label "wf_teloseq"
    cpus 1
    memory 2.GB
    input:
        tuple val(meta), path("staging_dir/*"), val(dirname)
    output:
        tuple val(meta), path(dirname)
    script:
    """
    mv staging_dir $dirname
    """
}


// workflow module
workflow pipeline {
    take:
        samples
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        metadata = samples.map { meta, reads, stats -> meta }.toList()

        // Read file and create metadata tuple
        Path ref = file(params.reference ?: "$projectDir/test_data/HG002qpMP_reference.fasta.gz", checkIfExists: true)

        //enzyme cut site to use
        enzyme_cut_site = params.enzyme_cut

        // Pass the channel to the check_reference process
        //check_reference(ref_ch)
        check_reference(ref,enzyme_cut_site)
        
        //remove duplicate reads - should not be common but just in case
        dedup = rmdup(samples.map{ meta, reads, stats -> [ meta,reads ] })

        //filter for telomere containing reads for 10 repeats within 60-500 bp of read
        telomeres = filter_telomeres(dedup)

        if (params.doublestranded) {
            //if both strands
            reversedreads=reversecomplement(dedup)
            telomeres2 = filter_telomeres2(reversedreads)
            telomeres3=combinefastq(telomeres,telomeres2)
            telomeres4=rmdup2(telomeres3)
            //filter for no telomere sequence 5 repeats for the last 60bp. Limited to 60 as cutsite for 2 chr arms is only 80bp from telomere. 
            //This is to filter out telomere only reads with no subtelomere to use for mapping and avoid telomere fragments.
            nontelomeres = filter_nontelomeres(telomeres4)
            remove_short_telomeres = remove_short2(telomeres4)
        } else {
            nontelomeres = filter_nontelomeres(telomeres)
            remove_short_telomeres = remove_short2(telomeres)
        }

        //remove short reads
        remove_short_nontelomeres = remove_short1(nontelomeres)

        //get subtelomere length information for plot
        sub1 = subtelomere(remove_short_nontelomeres)

        //This identifies reads with basecalling error to remove from the pipeline
        filtered = filter_motifs(remove_short_telomeres)

        //filtered telomere read fastq and filtered telomere-subtelomere fastq
        t1 = filter_motifs_reads1(remove_short_telomeres.join(filtered))
        t2 = filter_motifs_reads2(remove_short_nontelomeres.join(filtered))

        trim_adapters(t2)

        //telomere lengths plot for raw filtered error reads
        raw_telomere_analysis(trim_adapters.out.fastqtrimmed)

        //get min coverage
        read_count = trim_adapters.out.fastqtrimmed.countLines()

        //hard coded minimum read count for clustering, I did calculate on 20% of chr arm but when cov low then minimum read number would lead to inflated contigs. Default is 8 reads
        //I worry if go too low then snps/indels lead to novel contigs but are the same chr arm.  3000 telomere reads should be ~32 per chr arm so 8 is 25% but this is number
        //of clustered reads so not directly relatable and its how well I separated out the reads which won't be 100, it might be 50%.
        cov = params.cov_4cluster

        //coverage of final chr arm needs to be at least 20% average coverage otherwise likely duplicate contig.
        cov_filter=coverage_calc(trim_adapters.out.fastqtrimmed)
        
        //get stats on raw reads, telomere containing reads, exclude telomere only reads. Put output into variable for report.
        read_stats1 = fastq_stats(dedup)
        read_stats2 = fastq_stats2(t1)
        read_stats3 = fastq_stats3(trim_adapters.out.fastqtrimmed)


        //output results to channel for copying
        ch_to_publish = Channel.empty()
        | mix(
            software_versions | map { [it, null] },
            workflow_params | map { [it, null] },
        )

        //add to output channel telomere reads
        ch_to_publish = ch_to_publish
        | mix(
            trim_adapters.out.fastqtrimmed 
            | map { meta, reads -> [reads, "${meta.alias}/reads"] }
            | transpose
        )

        //add to output channel raw telomere results csv
        ch_to_publish = ch_to_publish 
        | mix(
            raw_telomere_analysis.out.plotraw 
            | map { meta, csv -> [csv, "${meta.alias}/results"] }
            | transpose
        )

        //add to output channel raw telomere pdf plot
        ch_to_publish = ch_to_publish 
        | mix(
            raw_telomere_analysis.out.pdfraw 
            | map { meta, pdf -> [pdf, "${meta.alias}/plots"] }
            | transpose
        )

        // NON MAPPING ROUTE REPORT  
        if (params.skipmapping) {
            // get all the per sample results together
            ch_per_sample_results = samples
            | map { meta, reads, stats_dir -> [meta, stats_dir] }
            | join(read_stats1)
            | join(read_stats2)
            | join(read_stats3)
            | join(raw_telomere_analysis.out.plotraw)
            | join(raw_telomere_analysis.out.covraw)
            | join(sub1)

            // collect results into a directory for the sample directory to avoid collisions
            ch_results_for_report = ch_per_sample_results
            | map {
                meta = it[0]
                rest = it[1..-1]
                [meta, rest, meta.alias]
            }
            | collectFilesInDir
            | map { meta, dirname -> dirname }

            //make report html with all information
            mappingreport=false
            report = makeReport(
                metadata,  
                software_versions,
                workflow_params,
                ch_results_for_report | collect,
                mappingreport
            )
        } else {
            //MAPPING ARM PIPELINE

            //denovo reference route
            if (params.denovo) {
                //remove subtelomere for raw and reference and map, output read subsets based upon mapping.
                clusterset = mapAndSeparateR1(trim_adapters.out.fastqtrimmed, check_reference.out.ref1, cov)
                //separate out each subset read file to a tuple in the new channel
                split_clusters = clusterset
                .flatMap { tuple ->
                // //tuple[0] is the list of four values and tuple[1] is the list of paths
                metadata = tuple[0]
                paths = tuple[1]
                //Create a new tuple for each path with the same key
                paths.collect { path -> [metadata, path] }
                }
                //combine to a channel so the fastq file with each text file so runs each rather than just once with one fastq file
                clusterchannel = split_clusters.combine(trim_adapters.out.fastqtrimmed, by: 0)
                //group tuple should collect all results with key. This is the vsearch clustering step for each set of reads identified as chr arm
                clusterout = clusterAndExtractR1(clusterchannel, cov).groupTuple()    
                //collect by groupTuple the files from each meta key to collate the contigs into one reference
                combineRefR1(trim_adapters.out.fastqtrimmed.join(clusterout), cov)

                /////////////////////////////////////////////////////////////////
                //phase 1 of clustering and assembling denovo reference done
                ////////////////////////////////////////////////////////////////


                //second round of phasing reads to cluster
                phaseagain = mapAndSeparateR2(combineRefR1.out.bam, cov)
                //separate out each file to a tuple in the new channel
                split_clusters2 = phaseagain
                .flatMap { tuple ->
                // //tuple[0] is the list of four values and tuple[1] is the list of paths
                metadata = tuple[0]
                paths = tuple[1]
                //Create a new tuple for each path with the same key
                paths.collect { path -> [metadata, path] }
                }
                //combine to a channel so the fastq file with each text file so runs each rather than just once with one fastq file
                clusterchannel2 = split_clusters2.combine(trim_adapters.out.fastqtrimmed, by: 0)
                //group tuple should collect all results with key. This is the vsearch clustering step for each set of reads identified as chr arm
                clusterout2 = clusterAndExtractR2(clusterchannel2, cov).groupTuple()  
                //collect by groupTuple the files from each meta key to collate the contigs into one reference
                combineRefR2(
                    trim_adapters.out.fastqtrimmed
                    | join(clusterout2)
                    | join(combineRefR1.out.ref1)
                    | join(cov_filter),
                    check_reference.out.ref1,
                    enzyme_cut_site,
                )

                /////////////////////////////////////////////////////////////////
                //phase 2 of clustering and assembling denovo reference done
                ////////////////////////////////////////////////////////////////


                //second round of phasing reads to cluster
                phaseagain2 = mapAndSeparateR3(combineRefR2.out.bam, cov)
                //separate out each file to a tuple in the new channel
                split_clusters3 = phaseagain2
                .flatMap { tuple ->
                // //tuple[0] is the list of four values and tuple[1] is the list of paths
                metadata = tuple[0]
                paths = tuple[1]
                //Create a new tuple for each path with the same key
                paths.collect { path -> [metadata, path] }
                }


                //combine to a channel so the fastq file with each text file so runs each rather than just once with one fastq file
                clusterchannel3 = split_clusters3.combine(trim_adapters.out.fastqtrimmed, by: 0)

                //group tuple should collect all results with key. This is the vsearch clustering step for each set of reads identified as chr arm
                clusterout3 = clusterAndExtractR3(clusterchannel3, cov).groupTuple()  

                //collect by groupTuple the files from each meta key to collate the contigs into one reference
                combineRefR3(
                    trim_adapters.out.fastqtrimmed
                    | join(clusterout3)
                    | join(combineRefR2.out.ref1)
                    | join(cov_filter),
                    check_reference.out.ref1,
                    enzyme_cut_site,
                )

                /////////////////////////////////////////////////////////////////
                //phase 3 of clustering and assembling denovo reference done
                ////////////////////////////////////////////////////////////////


                //ref1_paths = combineRefR2.out.ref1.map { path -> [meta, path] }     
                //ref1_paths = combineRefR3.out.ref1.map { tuple -> tuple[1] }
                //map filtered telomere reads to genome and filter using mapq (default=10)
                mappingbam(trim_adapters.out.fastqtrimmed, combineRefR3.out.ref1)
                //last telomere repeat location on the reference from the enzyme for each chr
                telomere_sites(combineRefR3.out.ref1)
                //filter bam with high, low and no stringency but including mapping quality filter applied in previous step
                filtering(mappingbam.output.alignments, combineRefR3.out.cutbed, telomere_sites.out.telomerebed)
                //get final telomere stats
                results(
                    filtering.out.combined
                    | join(raw_telomere_analysis.out.covraw)
                    | join(cov_filter),
                    combineRefR3.out.ref1,
                )

                //add to output channel
                ch_to_publish = ch_to_publish 
                | mix(
                    combineRefR3.out.ref1 
                    | map { meta, ref1  -> [[ref1], "$meta.alias/reference"] }
                    | transpose
                )
            } else {
                //using reference provided rather than de novo route.
                if (params.curation) {
                    //merge and correct de novo reference and new contigs
                    contigstoadd=file(params.curatedContigs, checkIfExists: true)
                    denovoref1=file(params.denovoRef, checkIfExists: true)
                    // TODO: do we need to make sure that `params.curatedContigs` and
                    // `params.denovoRef` are present when `params.curation`? if so, we
                    // should do this in the schema
                    manualCuration(
                        trim_adapters.out.fastqtrimmed
                        | join(cov_filter),
                        contigstoadd,
                        denovoref1,
                        check_reference.out.ref1,
                        enzyme_cut_site,
                    )
                    //last telomere repeat location on the reference from the enzyme for each chr
                    //ref1_paths = manualCuration.out.ref1.map { tuple -> tuple[1] }
                    telomere_sites(manualCuration.out.ref1)
                    //first cutsite location on the reference from the enzyme for each chr
                    cut_sites(manualCuration.out.ref1,enzyme_cut_site)
                    //map filtered telomere reads to genome and filter using mapq (default=10)
                    mappingbam(trim_adapters.out.fastqtrimmed, manualCuration.out.ref1)
                    //filter bam with high, low and no stringency but including mapping quality filter applied in previous step
                    filtering(mappingbam.output.alignments, cut_sites.out.cutbed, telomere_sites.out.telomerebed)
                    //get final telomere stats
                    results(
                        filtering.out.combined
                        | join(raw_telomere_analysis.out.covraw)
                        | join(cov_filter),
                        manualCuration.out.ref1
                    )

                    //add to output channel
                    ch_to_publish = ch_to_publish 
                        | mix(
                        manualCuration.out.ref1 
                        | map { meta, ref1  -> [[ref1], "$meta.alias/reference"] }
                        | transpose
                    )
                } else {
                    //last telomere repeat location on the reference from the enzyme for each chr
                    telomere_sites(check_reference.out.ref1)
                    //first cutsite location on the reference from the enzyme for each chr
                    cut_sites(check_reference.out.ref1,enzyme_cut_site)
                    //map filtered telomere reads to genome and filter using mapq (default=10)
                    mappingbam(trim_adapters.out.fastqtrimmed, check_reference.out.ref1)
                    //filter bam with high, low and no stringency but including mapping quality filter applied in previous step
                    filtering(mappingbam.output.alignments, cut_sites.out.cutbed, telomere_sites.out.telomerebed)
                    //get final telomere stats
                    results(
                        filtering.out.combined
                        | join(raw_telomere_analysis.out.covraw)
                        | join(cov_filter),
                        check_reference.out.ref1,
                    )
                }
            }

            // get all the per sample results together
            ch_per_sample_results = samples
            | map { meta, reads, stats_dir -> [meta, stats_dir] }
            | join(read_stats1)
            | join(read_stats2)
            | join(read_stats3)
            | join(raw_telomere_analysis.out.plotraw) 
            | join(raw_telomere_analysis.out.covraw)
            | join(sub1)
            | join(results.out.for_report)

            // collect results into a directory for the sample directory to avoid collisions
            ch_results_for_report = ch_per_sample_results
            | map {
                meta = it[0]
                rest = it[1..-1]
                [meta, rest, meta.alias]
            }
            | collectFilesInDir
            | map { meta, dirname -> dirname }

            //make report html file with all information
            mappingreport=true
            report = makeReport(
                    metadata,
                    software_versions,
                    workflow_params,
                    ch_results_for_report | collect,
                    mappingreport
                )

            //add to output channel, bam alignment files
            ch_to_publish = ch_to_publish
                | mix(
                    mappingbam.out.alignments
                    | map { meta, bam, bai -> [[bam, bai], "$meta.alias/alignments"] }
                    | transpose
                )

            //add to output channel, mapped csv files
            ch_to_publish = ch_to_publish 
                | mix(
                results.out.alldata 
                | map { meta, csv1,csv2,csv3,csv4,csv5 ,csv6 -> [[csv1, csv2, csv3, csv4, csv5, csv6], "${meta.alias}/results"] }
                | transpose
            )

            //add to output channel, mapped pdf plot files
            ch_to_publish = ch_to_publish 
                | mix(
                results.out.pdf 
                | map { meta, pdf1,pdf2,pdf3,pdf4,pdf5 ,pdf6 -> [[pdf1, pdf2, pdf3, pdf4, pdf5, pdf6], "${meta.alias}/plots"] }
                | transpose
            )
            //add to output channel, reference used for mapping final results
            ch_to_publish = ch_to_publish 
                | mix(
                mappingbam.out.mappingref 
                | map { meta, mappingref  -> [[mappingref], "$meta.alias/alignments"] }
                | transpose
            )
            
            //add to output channel, filtered strict bam files
            ch_to_publish = ch_to_publish 
                | mix(
                filtering.out.finalbam 
                | map { meta, bam, bai  -> [[bam, bai], "$meta.alias/alignments"] }
                | transpose
            )

            //add to output channel, filtered lenient bam files
            ch_to_publish = ch_to_publish 
                | mix(
                filtering.out.lowfinalbam 
                | map { meta, bam, bai  -> [[bam, bai], "$meta.alias/alignments"] }
                | transpose
            )
        }

    //this emits the report, files to output directory and telemetry information
    emit:
        report
        combined_results_to_publish = ch_to_publish
        workflow_params
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    // demo mutateParam
    if (params.containsKey("mutate_fastq")) {
        CWUtil.mutateParam(params, "fastq", params.mutate_fastq)
    }

    Map ingress_args = [
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "stats": true,   // TODO: we might wanna use these instead of the seqkit stats
        "fastcat_extra_args": "",   // TODO: we could use this to filter based on read length + quality
    ]
    if (params.fastq) {
        samples = fastq_ingress(ingress_args + [
            "input":params.fastq
        ])
    } else {
        // if we didn't get a `--fastq`, there must have been a `--bam` (as is codified
        // by the schema)
        samples = xam_ingress(ingress_args + [
            "input":params.bam,
            "return_fastq": true,
            "keep_unaligned": true,
        ])
    }

    pipeline(samples)

    // publish results
    pipeline.out.combined_results_to_publish
    | toList
    | flatMap | concat (
        pipeline.out.report.concat(pipeline.out.workflow_params)
        | map { [it, null] }
    )
    | output

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
