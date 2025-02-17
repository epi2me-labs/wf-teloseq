#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { mapping } from './modules/local/common'
include { fastq_stats as fastq_stats_all } from './modules/local/common'
include { fastq_stats as fastq_stats_telo } from './modules/local/common'
include { fastq_stats as fastq_stats_telo_subtelo } from './modules/local/common'
include { process_reads } from './modules/local/process'
include { check_reference; trim_restriction_site; trim_restriction_site_ref; coverage_filter2; nm_filter; telomere_sites; filtering; results } from './modules/local/to_review'

process getVersions {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    output: path "versions.txt"
    script:
    """
    python -c "import matplotlib as mpl; print(f'matplotlib,{mpl.__version__}')" >> versions.txt
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas as pd; print(f'pandas,{pd.__version__}')" >> versions.txt
    python -c "import numpy as np; print(f'numpy,{np.__version__}')" >> versions.txt
    python -c "import ezcharts; print(f'ezcharts,{ezcharts.__version__}')" >> versions.txt
    seqkit version | grep "eqkit" | awk '{print "seqkit," \$2}' >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    blastn -version | awk 'NR==1 {print "blastn," \$2; exit}' >> versions.txt
    csvtk version | awk '{print "csvtk," \$2}' >> versions.txt
    minimap2 --version | awk '{print "minimap2," \$1}' >> versions.txt
    cutadapt --version | awk '{print "cutadapt," \$1}' >> versions.txt
    ( seqtk 2>&1 || true ) | grep "Version" | awk '{print "seqtk," \$2}' >> versions.txt
    """
}

// Very unsure about this fixed 5
process coverage_calc {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    input:
        tuple val(meta), path("telomere.fastq")
    output:
        tuple val(meta), stdout, emit: cov
    script:
    """
    if [[ "${params.min_coverage}" != -1 ]]; then
        coverage=${params.min_coverage}
    else
        read_count=\$(wc -l < telomere.fastq | awk '{print \$1}')
        coverage=\$(( read_count / 4 / ${params.exp_chr_num} * ${params.min_coverage_percent} / 100 ))
        if (( coverage < 5 )); then
            coverage=5
        fi
    fi
    echo -n \$coverage
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
        // Define default reference at workflow level, replaced by user option if supplied.
        def default_ref = params.reference ? file(params.reference, checkIfExists: true) : file("$projectDir/data/HG002qpMP_reference.fasta.gz", checkIfExists: true)
        
        // Get software and workflow paramaters for later report
        workflow_params = getParams()
        software_versions = getVersions()

        // Filter reads down to high quality telomere repeat spanning reads. Errors if no reads left after processing all samples
        // Saves some length based statistics about sample read datasets
        filtered_samples = process_reads(samples)

        valid_samples_with_stats = filtered_samples
            .filter { _meta, sub_telomeric_fastq_file, _length_sum, _length_records, _stats_dir ->
                sub_telomeric_fastq_file.size() > 0
            }.map { meta, _sub_telomeric_fastq_file, _length_sum, _length_records, _stats_dir  ->
                if (!meta.containsKey('reference')) {
                    meta.reference = default_ref
                }
                [meta, _sub_telomeric_fastq_file, _length_sum, _length_records, _stats_dir]
            }
            .ifEmpty {
                throw new Exception("No valid samples found. Exiting workflow.")
            }

    emit:
        metadata = valid_samples_with_stats.map { items -> items[0] }.collect()
        software_versions
        workflow_params
        valid_samples_with_stats
}

/**
 * Workflow: prealignment_stats
 * 
 * Calculates some final statistics for valid sequencing samples, 
 * including expected coverage and seqkit stats statistics for subtelomere containing reads.
 *
 */
workflow prealignment_stats {
    take:
    valid_samples_with_stats
    main:
        // Calculate the expected coverage
        // Calculate the expected coverage of telomere reads per chr arm at 15% of average but can be overridden by min_coverage
        meta_and_reads = valid_samples_with_stats.map{items -> {
            def meta = items[0]
            def reads = items[1]
            [meta, reads]
        }}
        cov_filter = coverage_calc(meta_and_reads)
        
        // Statistics gathering
        // Read seqkit stats on trimmed telomere and subtelomere containing reads  
        subtelomere_stats = fastq_stats_telo_subtelo(meta_and_reads, "subtelomere_reads_stats.txt", "Telomere_subtelomere_Reads")
        
        // Prepare final results channel, to be added to later depending upon path, removing read file path
        per_sample_stats = valid_samples_with_stats
            | map { meta, _reads, length_sum, length_records, stats_dir -> [meta, length_sum, length_records, stats_dir] }
            | join(subtelomere_stats)


    emit:
        per_sample_stats
        // Couldn't the join above allow us to access this file anyway?
        cov_filter
    
}

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
        valid_samples_with_stats = preprocessing.valid_samples_with_stats

        pre_stats = prealignment_stats(valid_samples_with_stats)
        per_sample_stats = pre_stats.per_sample_stats
        cov_filter = pre_stats.cov_filter

        // Link sample to a reference
        reference_channel = valid_samples_with_stats
            .map { items ->
                def meta = items[0]
                def ref = file(meta.reference)
                return tuple(meta, ref)
            }


        if (params.skip_mapping) {
            // Collect results into directories to avoid collisions
            ch_results_for_report = per_sample_stats
                | map { items ->
                    def meta = items[0]
                    def rest = items[1..-1]
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
            combined_results_to_publish = Channel.empty()
            | mix(
                software_versions | map { it -> [it, null] },
                workflow_params | map {it -> [it, null] }
            )
            | mix(
                // Telomere reads
                valid_samples_with_stats.map { items -> 
                    def reads = items[1]
                    def meta = items[0]
                    [reads, "${meta.alias}/reads"] }
                .transpose(),
                // Raw telomere results CSV
                per_sample_stats.map { meta, _length_summary, lengths_records, _stats_dir, _seqkit_stats -> [lengths_records, "${meta.alias}/results"] }.transpose()
            )

        } else {
            ////////////////////////////////////////////
            // MAPPING ARM PIPELINE
            ////////////////////////////////////////////

            // Prepare reference by identifying telomere contigs and extracting/orientating
            check_reference(reference_channel)

            // Trim reads by restriction site
            processed_read_channel = trim_restriction_site(valid_samples_with_stats.map{items -> items[0..1]})

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
                | join(per_sample_stats.map{items -> [items[0], items[1]]}, by: 0)
                | join(cov_filter, by: 0)
                | join(check_reference.out.reference, by: 0)
            )

            // Collect all per-sample results
            ch_per_sample_results = per_sample_stats
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
            combined_results_to_publish = Channel.empty()
                | mix(
                    software_versions | map { it -> [it, null] },
                    workflow_params | map {it -> [it, null] }
                )
                | mix(
                    // Telomere reads
                    valid_samples_with_stats.map { items -> 
                        def reads = items[1]
                        def meta = items[0]
                        [reads, "${meta.alias}/reads"] }
                    .transpose(),
                    // Raw telomere results CSV
                    per_sample_stats.map { meta, _length_summary, lengths_records, _stats_dir, _seqkit_stats -> [lengths_records, "${meta.alias}/results"] }.transpose()
                )

            if (!params.skip_mapping) {

                combined_results_to_publish = combined_results_to_publish | mix(
                    // Mapped CSV files
                    results.out.for_report.map { items -> 
                        def meta = items[0]
                        def files = items[1..-1]
                        [files, "${meta.alias}/results"]
                    }.transpose(),
                    // Reference used for mapping
                    coverage_filter2.out.reference.map { meta, reference -> [[reference], "${meta.alias}/alignment"] }.transpose(),
                    // Filtered tagged BAM files
                    filtering.out.finalbam.map { meta, bam, bai -> [[bam, bai], "${meta.alias}/alignment"] }.transpose()
                )
            }
        }


    emit:
        report
        combined_results_to_publish
        telemetry = workflow_params
        workflow_params

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
        | map { it -> [it, null] }
    )
    | publish


}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
