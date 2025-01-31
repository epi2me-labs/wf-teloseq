#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { mapping } from './modules/local/common'
include { fastq_stats as fastq_stats_all } from './modules/local/common'
include { fastq_stats as fastq_stats_telo } from './modules/local/common'
include { fastq_stats as fastq_stats_telo_subtelo } from './modules/local/common'
include { filter_motifs_reads as filter_motifs_reads1 } from './modules/local/common'
include { filter_motifs_reads as filter_motifs_reads2 } from './modules/local/common'


OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process getVersions {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    output: path "versions.txt"
    script:
    """
    python -c "import matplotlib as mpl; print(f'matplotlib,{mpl.__version__}')" >> versions.txt
    python -c "import pyfastx; print(f'pyfastx,{pyfastx.__version__}')" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas as pd; print(f'pandas,{pd.__version__}')" >> versions.txt
    python -c "import numpy as np; print(f'numpy,{np.__version__}')" >> versions.txt
    python -c "import ezcharts; print(f'ezcharts,{ezcharts.__version__}')" >> versions.txt
    seqkit version | grep "eqkit" | awk '{print "seqkit," \$2}' >> versions.txt
    python -c 'from importlib.metadata import version; v=version("edlib"); print(f"edlib,{v}")' \
        >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    blastn -version | awk 'NR==1 {print "blastn," \$2; exit}' >> versions.txt
    csvtk version | awk '{print "csvtk," \$2}' >> versions.txt
    minimap2 --version | awk '{print "minimap2," \$1}' >> versions.txt
    cutadapt --version | awk '{print "cutadapt," \$1}' >> versions.txt
    vsearch --version 2>&1 | grep "vsearch " | sed 's/,.*//' | sed 's/ /,/' | sed 's/_.*//' >> versions.txt
    ( seqtk 2>&1 || true ) | grep "Version" | awk '{print "seqtk," \$2}' >> versions.txt
    """
}

// Note: First open and operation on input reads
process checkReads {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq"), path(stats)
    output:
        tuple val(meta), path("reads.fastq"), env(IS_VALID)
    script:
    """
    read_count=\$(seqkit grep -s -R 60:500 -P -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" reads.fastq | seqkit stats -a | awk 'NR==2 {print \$4}' | tr -d ',')
    
    if [[ \$read_count -gt 0 ]]; then
        IS_VALID=true
    else
        IS_VALID=false
    fi
    """
}
// Looks for telomeric repeats in the fastq and filters those reads into a file
// Note: Second open and operation on input reads
process filter_reads {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere.fastq")
    script:
    """
    seqkit grep -s -R 60:500 -P -p "TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC" reads.fastq > telomere.fastq
    """
}
// Looks for telomeric repeats in the last 70bp. These reads are considered short, as the cut site is after this.
// These reads are added into their own file
// Note: Third open and operation on input reads
process filter_nontelomeres {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("telomere_and_non_telomere.fastq")
    // filter further subset that should not have telomere sequence in the last 70bp as shortest cut site is beyond this so should have atleast this amount of sequence at end that
    // is not telomere and don't want telomere only reads as could be fragments or artefacts. 
    script:
    """
    seqkit grep -v -s -R -70:-1 -P -p "TAACCCTAACCCTAACCCTAACCC" reads.fastq > telomere_and_non_telomere.fastq
    """
}
// Note: Third open and operation on input reads
// Cuts a third different length of the telomere repeat from the 5' end
// Then reads the ouput fastq and writes a file of read lengths
process subtelomere {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("subtelomere.txt")
    script:
    """
    cutadapt -g "TAACCCTAACCCTAACCCTAACCCTAACCC;rightmost" -e 0 -o subtelomere.fastq reads.fastq
    awk 'NR%4 == 2 {print length(\$0)}' subtelomere.fastq | tr -s ' ' '\t' > subtelomere.txt
    """
}

// Note: fourth open of input reads
// AFAIK - this is finding the location of specifically common basecalling error motifs, and also the location of the last telomere repeat in the read
// It is then using AWK to get the last telomere location, and passing both into a python script that dumps those read ids into a file if 
// there are a certain number of basecalling errors within a window, AFTER the last telomere repeat.
process filter_motifs {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("remove_read_ids.txt")
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
        last_telomere_locations.tsv \
        error_locations.tsv \
        remove_read_ids.txt \
        $params.filter_error_motifs_window_size \
        $params.filter_error_motifs_max_count
    """
}

process trim_adapters {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), path("reads_trimmed.fastq"), emit: fastqtrimmed
    //trim adapter. why not use cutadapt, because adapter not basecalled well and precision important to end motifs
    script:
    """
    workflow-glue trim_adapters reads.fastq $params.adapters_to_trim > reads_trimmed.fastq
    """
}


process telomere_lengths_raw {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        tuple val(meta), path("telomere.fastq")
    output:
        tuple val(meta), path("sample_raw_coverage.csv"), emit: covraw
        tuple val(meta), path("sample_raw_per_read_telomere_length.csv"), emit: plotraw
    script:
    // Note: this assumes `telomere length == coordinate of telomere boundary` (i.e.
    // that there is no sequence that is not telomere at the beginnings of the reads)
    """
    workflow-glue get_telomere_boundaries_in_fastx telomere.fastq \\
        > telomere_read_length.txt
    workflow-glue summarise_telomere_lengths --mode raw --telomere_lengths telomere_read_length.txt --output_prefix sample
    """
}

process coverage_calc {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("reads.fastq")
    output:
        tuple val(meta), stdout, emit: cov
    script:
    """
    if [[ "${params.min_coverage}" != -1 ]]; then
        coverage=${params.min_coverage}
    else
        read_count=\$(wc -l < reads.fastq | awk '{print \$1}')
        coverage=\$(( read_count / 4 / ${params.exp_chr_num} * ${params.min_coverage_percent} / 100 ))
        if (( coverage < 5 )); then
            coverage=5
        fi
    fi
    echo -n \$coverage
    """
}


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


// Make report html file using staged file dir and other files
process makeReport {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    input:
        val meta_array
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
        def simplified_meta = meta_array.collect { meta ->
            [
                alias: meta.alias,
                barcode: meta.barcode,
                type: meta.type,
                run_ids: meta.run_ids,
                basecall_models: meta.basecall_models,
                is_unaligned: meta.is_unaligned,
                n_seqs: meta.n_seqs
            ]
        }
        def json_str = new groovy.json.JsonBuilder(simplified_meta).toPrettyString()
        """
        echo '${json_str}' > metadata.json

        workflow-glue report $report_name \\
            --metadata metadata.json \\
            --versions versions \\
            --params params.json \\
            --data data \\
            $extra_arg
        """
}

// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process publish {
    // publish inputs to output directory
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    script:
        """
        echo "Writing output files"
        """
}

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


workflow generic_preprocessing {
    take:
        samples
    
    main:
        // Get software and workflow paramaters for later report
        software_versions = getVersions()
        workflow_params = getParams()
        
        // Check samples for validity, otherwise a sample without telomere reads could crash pipeline
        checked_samples = checkReads(samples)
        
        // Define default reference at workflow level, replaced by user option if supplied.
        def default_ref = params.reference ? file(params.reference, checkIfExists: true) : file("$projectDir/data/HG002qpMP_reference.fasta.gz", checkIfExists: true)
        
        // Process valid samples and add reference key to meta data if not supplied from sample_sheet parameter
        valid_samples = checked_samples
            .filter { meta, reads, is_valid -> 
                def keep = is_valid.toBoolean()
                return keep
            }
            .map { meta, reads, is_valid ->
                if (!meta.containsKey('reference')) {
                    meta.reference = default_ref
                }
                return tuple(meta, reads)
            }

        // Create and verify reference channel using sample meta data linking to reference
        reference_channel = valid_samples
            .map { meta, reads ->
                def ref = file(meta.reference)
                return tuple(meta, ref)
            }
        
        // Process stats
        stats_channel = samples
            .map { meta, reads, stats_dir -> 
                tuple(meta.alias, stats_dir) 
            }
        
        // Join samples with stats for later report
        valid_samples_with_stats = valid_samples
            .map { meta, reads -> 
                tuple(meta.alias, meta, reads)
            }
            .combine(stats_channel, by: 0)
            .map { alias, meta, reads, stats_dir -> 
                tuple(meta, reads, stats_dir) 
            }
            
        // Initial filtering and processing
        // Identifying telomere containing reads
        telomere_set = filter_reads(valid_samples)
        // Identifying telomere containing reads that also are not just telomere but have subtelomere    
        nontelomeres = filter_nontelomeres(telomere_set)
        // For each reads calculating subtelomere length for later report plot and sample debug    
        sub1 = subtelomere(nontelomeres)
        // Filter out reads with sequencing error clusters of x length    
        filtered = filter_motifs(nontelomeres)
        // Filtered telomere containing reads
        t1 = filter_motifs_reads1(telomere_set.join(filtered))
        // Filtered telomere reads that contain subtelomere
        t2 = filter_motifs_reads2(nontelomeres.join(filtered))
        // Trim the adapters off the reads
        trimmed_fastq_channel = trim_adapters(t2)
        // Get the raw per read telomere lengths before mapping    
        (telomere_lengths_raw_channel, telomere_lengths_plotraw_channel) = telomere_lengths_raw(trimmed_fastq_channel)
        // Calculate the expected coverage of telomere reads per chr arm at 15% of average but can be overridden by min_coverage
        cov_filter = coverage_calc(trimmed_fastq_channel)
        
        // Statistics gathering
        // Read seqkit stats on raw fastq files
        read_stats1 = fastq_stats_all(valid_samples, "raw.txt", "Raw_Reads" )
        // Read seqkit stats on telomere identified reads  
        read_stats2 = fastq_stats_telo(t1, "telomere.txt", "Telomere_Reads")
        // Read seqkit stats on trimmed telomere and subtelomere containing reads  
        read_stats3 = fastq_stats_telo_subtelo(trimmed_fastq_channel, "telomere_subtelomere.txt", "Telomere_subtelomere_Reads")
        
        // Prepare final results channel, to be added to later depending upon path
        ch_per_sample_results = valid_samples_with_stats
            | map { meta, reads, stats_dir -> [meta, stats_dir] }
            | join(read_stats1)
            | join(read_stats2)
            | join(read_stats3)
            | join(telomere_lengths_raw.out.plotraw)
            | join(telomere_lengths_raw.out.covraw)
            | join(sub1)
            
    emit:
        valid_samples_with_stats
        metadata = valid_samples.map { meta, reads -> meta }.collect()
        software_versions
        workflow_params
        trimmed_fastq_channel
        ch_per_sample_results
        cov_filter
        telomere_lengths_raw_channel
        telomere_lengths_plotraw_channel
        reference_channel
}


// workflow module
workflow pipeline {
    take:
        samples

    main:
        // Invoke the generic_preprocessing workflow
        preprocessing = generic_preprocessing(samples)
        
        // Extract outputs from preprocessing
        software_versions = preprocessing.software_versions
        workflow_params = preprocessing.workflow_params
        metadata = preprocessing.metadata
        trimmed_fastq_channel = preprocessing.trimmed_fastq_channel
        ch_per_sample_results = preprocessing.ch_per_sample_results
        cov_filter = preprocessing.cov_filter
        telomere_lengths_raw_channel = preprocessing.telomere_lengths_raw_channel
        telomere_lengths_plotraw_channel = preprocessing.telomere_lengths_plotraw_channel
        reference_channel = preprocessing.reference_channel
        
        if (params.skip_mapping) {
            // Collect results into directories to avoid collisions
            ch_results_for_report = ch_per_sample_results
                | map {
                    def meta = it[0]
                    def rest = it[1..-1]
                    [meta, rest, meta.alias]
                }
                | collectFilesInDir
                | map { _meta, dirname -> dirname }

            // Generate the report
            mappingreport=false
            report = makeReport(
                metadata,  
                software_versions,
                workflow_params,
                ch_results_for_report | collect,
                mappingreport
            )

        // Initialize ch_to_publish as an empty channel since mapping is skipped
        ch_to_publish = Channel.empty()

        } else {
            ////////////////////////////////////////////
            // MAPPING ARM PIPELINE
            ////////////////////////////////////////////

            // Prepare reference by identifying telomere contigs and extracting/orientating
            check_reference(reference_channel)

            // Trim reads by restriction site
            processed_read_channel = trim_restriction_site(trimmed_fastq_channel)

            ////////////////////////////////////////////
            // Map to Reference and Generate Results
            ////////////////////////////////////////////

            reference_channel2 = check_reference.out.reference

            // Trim reads by restriction site for reference
            trim_restriction_site_ref(reference_channel2)

            //last telomere repeat location on the reference from the enzyme for each chr
            telo_site_channel = telomere_sites(trim_restriction_site_ref.out.reference)

            //map filtered telomere reads to genome and filter using mapq (default=4)
            mapping(processed_read_channel.join(trim_restriction_site_ref.out.reference))

            //Remove incorrectly mapped reads via NM
            nm_filter(mapping.out.alignment)

            //filter for min sequences 10 (0.15 of average coverage).
            coverage_filter2(nm_filter.out.alignment
            | join (trim_restriction_site_ref.out.reference)
            | join (cov_filter))

            //filter bam with high, low and no stringency but including mapping quality filter applied in previous step
            filtering(coverage_filter2.out.alignment
            | join (trim_restriction_site_ref.out.cutsites)
            | join (telo_site_channel),
            params.telomere_margin,
            )

            // Generate final telomere statistics
            results(
                filtering.out.finalbam
                | join(telomere_lengths_raw_channel, by: 0)
                | join(cov_filter, by: 0)
                | join(check_reference.out.reference, by: 0)
            )

            // Collect all per-sample results
            ch_per_sample_results = ch_per_sample_results
                | join(results.out.for_report, by: 0)

            // Collect results into directories to avoid file collisions
            ch_results_for_report = ch_per_sample_results
                | map {
                    def meta = it[0]
                    def rest = it[1..-1]
                    def dirname = meta.alias  // Use the alias string for the directory name
                    [meta, rest, dirname]  // Keep meta object but use alias for dirname
                }
                | collectFilesInDir
                | map { meta, dirname -> dirname }

            // Generate the mapping version of report
            mappingreport=true
            report = makeReport(
                metadata,
                software_versions,
                workflow_params,
                ch_results_for_report | collect,
                mappingreport
            )


            ////////////////PUBLISH RESULTS SECTION///////////////////////////////

            // Initialize channel for publishing results
            ch_to_publish = Channel.empty()
                | mix(
                    software_versions | map { [it, null] },
                    workflow_params | map { [it, null] }
                )
                | mix(
                    // Telomere reads
                    trimmed_fastq_channel.map { meta, reads -> [reads, "${meta.alias}/reads"] }.transpose(),
                    // Raw telomere results CSV
                    telomere_lengths_plotraw_channel.map { meta, csv -> [csv, "${meta.alias}/results"] }.transpose()
                )

            if (!params.skip_mapping) {

                ch_to_publish = ch_to_publish | mix(
                    // Mapped CSV files
                    results.out.for_report.map {
                        meta = it[0]
                        files = it[1..-1]
                        [files, "${meta.alias}/results"]
                    }.transpose(),
                    // Reference used for mapping
                    coverage_filter2.out.reference.map { meta, reference -> [[reference], "${meta.alias}/alignment"] }.transpose(),
                    // Filtered tagged BAM files
                    filtering.out.finalbam.map { meta, bam, bai -> [[bam, bai], "${meta.alias}/alignment"] }.transpose()
                )
            }
        }

    // Emit outputs: report, results to publish, and telemetry
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
        "stats": true,
        "fastcat_extra_args": "--min_qscore ${params.read_quality} --min_length ${params.min_length}", 
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
    | publish

}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
