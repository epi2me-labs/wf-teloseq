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
        path("reference")  // can be compressed or uncompressed
        val(enzyme_motif)
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
    # TODO: what does the above comment mean?

    workflow-glue extract_reference \\
        reference \\
        TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC \\
        $enzyme_motif \\
        300 \\
        > extracted.fasta

    # TODO: the '300' above and '4000' below should be exposed as advanced param
    workflow-glue extend_telomere extracted.fasta TAACCC 4000 > reference.fasta
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
    workflow-glue trim_adapters reads.fastq $params.adapters_to_trim > reads_trimmed.fastq
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
    script:
    String telomere_pattern = "TAACCCTAACCCTAACCCTAACCCTAACCC"
    //seqkit forward strand the telomere repeat and identify all locations for eachr read
    //reverse the file so last telomere location for each read is first and print one row for each read, effectively getting last telomere match position for each read
    //retrieve just id and length information I need
    //search for basecalling error in telomere
    //Idnetify reads with high intensity so lots of kmer error clustered error within the telomere as this will not map and get softclipped, plus not useful for telomere length.
    //ensure just one read ID and @ is removed as not needed later on.
    """
    # `seqkit locate` output is grouped by input sequence; we can use `awk` to print the
    # last line of each group (i.e. read in our case)
    seqkit locate -j 1 -P -m 1 -p $telomere_pattern reads.fastq \
    | awk -F '\\t' '
        BEGIN { OFS="\\t" }
        {
            if (\$1!=prev_id && prev_id) print prev_id, prev_start
            prev_id=\$1
            prev_start=\$5
        }
        END { print \$1, \$5 }
    ' > last_telomere_locations.tsv

    seqkit locate -j 1 -P -m 0 -p $params.filter_error_motifs reads.fastq \
    | cut -f1,5 > error_locations.tsv

    workflow-glue error_reads \
        error_locations.tsv \
        last_telomere_locations.tsv \
        removereadids.txt
        # consider reinstating these parameters
        #--window_size $params.filter_error_motifs_window_size \
        #--max_count $params.filter_error_motifs_max_count
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
    # TODO: could do indexing during `samtools sort` to be a little faster
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
        tuple val(meta2), path("telomere_boundaries.bed"), emit: telomerebed
    script:
    """
    workflow-glue get_telomere_boundaries_in_fastx tel_reference.fasta \
        > telomere_boundaries.bed
    """
}

process filtering {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("align.bam"), path("align.bam.bai")
        tuple val(meta2), path("cutsites.bed")
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
    seqkit bam align.bam 2> sk-bam.tsv

    workflow-glue filter_bam_reads_output_id \\
        --seqkit_bam_out sk-bam.tsv \\
        --cut_sites_bed cutsites.bed \\
        --telomere_ends_bed telomeresites.txt \\
        --beyond_telomere_margin 80 \\
        --close_to_cutsite_margin 25 \\
        --output_strict_filter strict-keep-ids.txt \\
        --output_lenient_filter lenient-keep-ids.txt \\
        --output_no_filter no-extra-filter-keep-ids.txt

    samtools view -N strict-keep-ids.txt align.bam -o highfiltered.bam
    samtools index highfiltered.bam

    samtools view -N lenient-keep-ids.txt align.bam -o lowfiltered.bam
    samtools index lowfiltered.bam

    samtools view -N no-extra-filter-keep-ids.txt align.bam -o nofiltered.bam
    samtools index nofiltered.bam
    """
}

process telomere_lengths_raw {
    label "wf_teloseq"
    cpus   = 1
    memory 2.GB
    input:
        tuple val(meta), path("telomere.fastq")
    output:
        tuple val(meta), path("Sample_raw_Coverage.csv"), emit: covraw
        tuple val(meta), path("Sample_raw_Per_Read_telomere_length.csv"), emit: plotraw
    script:
    // Note: this assumes `telomere length == coordinate of telomere boundary` (i.e.
    // that there is no sequence that is not telomere at the beginnings of the reads)
    """
    workflow-glue get_telomere_boundaries_in_fastx telomere.fastq \\
        > telomere_read_length.txt
    workflow-glue summarise_raw_telomere_lengths telomere_read_length.txt Sample
    """
}

process results {
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
            path("lowfiltered_Per_Read_telomere_length.csv"),
            path("nofiltered_chr_arm_Coverage.csv"),
            path("nofiltered_Per_Read_telomere_length.csv"),
            path("output.csv"),
            emit: for_report
    script:
    //search for locations of telomere sequences (x5 repeats) in individual reads.
    //Reverse locations file to select last occurance of each telomere match, thereby selecting end position of telomere
    //remove the first location of telomere from the last to remove any adapter length contributing.
    """
    cov=\$(cat coverage.txt)
    samtools faidx mapping_ref.fasta

    seqkit bam highfiltered.bam -f Read,Ref,Acc,Strand,IsSec,IsSup \\
        2> highfiltered.sk-bam.tsv
    samtools fastq highfiltered.bam > highfiltered.fastq
    # TODO: could modify script to read either FASTQ or BAM
    # TODO: we should only call the script for finding the boundaries once on all the
    #       reads (and split them up based on the filter afterwards)
    workflow-glue get_telomere_boundaries_in_fastx highfiltered.fastq > high_telomere_read_length.txt
    workflow-glue summarise_mapped_telomere_lengths \\
        --seqkit-bam highfiltered.sk-bam.tsv \\
        --ref-fai mapping_ref.fasta.fai \\
        --min-coverage \$cov \\
        --telomere-lengths high_telomere_read_length.txt \\
        --output-prefix highfiltered

    seqkit bam lowfiltered.bam -f Read,Ref,Acc,Strand,IsSec,IsSup \\
        2> lowfiltered.sk-bam.tsv
    samtools fastq lowfiltered.bam > lowfiltered.fastq
    workflow-glue get_telomere_boundaries_in_fastx lowfiltered.fastq > low_telomere_read_length.txt
    workflow-glue summarise_mapped_telomere_lengths \\
        --seqkit-bam lowfiltered.sk-bam.tsv \\
        --ref-fai mapping_ref.fasta.fai \\
        --min-coverage \$cov \\
        --telomere-lengths low_telomere_read_length.txt \\
        --output-prefix lowfiltered

    seqkit bam nofiltered.bam -f Read,Ref,Acc,Strand,IsSec,IsSup \\
        2> nofiltered.sk-bam.tsv
    samtools fastq nofiltered.bam > nofiltered.fastq
    workflow-glue get_telomere_boundaries_in_fastx nofiltered.fastq > no_telomere_read_length.txt
    workflow-glue summarise_mapped_telomere_lengths \\
        --seqkit-bam nofiltered.sk-bam.tsv \\
        --ref-fai mapping_ref.fasta.fai \\
        --min-coverage \$cov \\
        --telomere-lengths no_telomere_read_length.txt \\
        --output-prefix nofiltered

    workflow-glue combine_result_CSVs \\
        raw_coverage.csv \\
        nofiltered.csv \\
        lowfiltered.csv \\
        highfiltered.csv \\
        output.csv
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
         // TODO: `nontelomeres` are telomeric reads that reach into the chromosomes
         // `telomeres` also contains reads that don't reach into the chromosomes (they
         // are only used for stats)
        }
        // TODO: we can use `bcftools consensus` now (wasn't used because of version
        // issues)

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
        telomere_lengths_raw(trim_adapters.out.fastqtrimmed)

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
            telomere_lengths_raw.out.plotraw 
            | map { meta, csv -> [csv, "${meta.alias}/results"] }
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
            | join(telomere_lengths_raw.out.plotraw)
            | join(telomere_lengths_raw.out.covraw)
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

            // Read file and create metadata tuple
            Path ref = file(params.reference ?: "$projectDir/data/HG002qpMP_reference.fasta.gz", checkIfExists: true)

            // prepare ref
            check_reference(ref, params.restriction_site)

            //last telomere repeat location on the reference from the enzyme for each chr
            telomere_sites(check_reference.out.ref1)
            //first cutsite location on the reference from the enzyme for each chr
            cut_sites(check_reference.out.ref1, params.restriction_site)
            //map filtered telomere reads to genome and filter using mapq (default=10)
            mappingbam(trim_adapters.out.fastqtrimmed, check_reference.out.ref1)
            //filter bam with high, low and no stringency but including mapping quality filter applied in previous step
            filtering(mappingbam.output.alignments, cut_sites.out.cutbed, telomere_sites.out.telomerebed)
            //get final telomere stats
            results(
                filtering.out.combined
                | join(telomere_lengths_raw.out.covraw)
                | join(cov_filter),
                check_reference.out.ref1,
            )

            // get all the per sample results together
            ch_per_sample_results = samples
            | map { meta, reads, stats_dir -> [meta, stats_dir] }
            | join(read_stats1)
            | join(read_stats2)
            | join(read_stats3)
            | join(telomere_lengths_raw.out.plotraw) 
            | join(telomere_lengths_raw.out.covraw)
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
                results.out.for_report
                | map {
                    meta = it[0]
                    files = it[1..-1]
                    [files, "${meta.alias}/results"]
                }
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
