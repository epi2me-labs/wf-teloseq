#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { align_and_process } from './modules/local/alignment'
include { fastq_stats; collectFilesInDir } from './modules/local/common'
include { process_reads } from './modules/local/process'
include { results } from './modules/local/to_review'

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

// Make report html file using staged file dir and other files
process makeReport {
    label "wf_teloseq"
    cpus   1
    memory '2 GB'
    publishDir (
        params.out_dir,
        mode: "move",
        overwrite: true
    )
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
            --workflow-version v${workflow.manifest.version} \\
            $extra_arg
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
        meta_and_reads = valid_samples_with_stats.map{items -> {
            def meta = items[0]
            def reads = items[1]
            [meta, reads]
        }}        
        // Statistics gathering
        // Read seqkit stats on trimmed telomere and subtelomere containing reads  
        subtelomere_stats = fastq_stats(meta_and_reads, "subtelomere_reads_stats.txt", "Telomere_subtelomere_Reads")
        
        // Prepare final results channel, to be added to later depending upon path, removing read file path
        per_sample_stats = valid_samples_with_stats
            | map { meta, _reads, length_sum, length_records, stats_dir -> [meta, length_sum, length_records, stats_dir] }
            | join(subtelomere_stats)


    emit:
        per_sample_stats
    
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

        // Link sample to a reference
        reference_channel = valid_samples_with_stats
            .map { items ->
                def meta = items[0]
                def ref = file(meta.reference)
                return tuple(meta, ref)
            }

        if (params.mapping) {
            //map filtered telomere reads to genome and filter using mapq (default=4)
            cleaned_alignments = align_and_process(valid_samples_with_stats.map{items -> 
                def meta = items[0]
                def reads = items[1]
                [meta, reads, file(meta.reference)]
            })

            // Generate final telomere statistics
            results(
                cleaned_alignments
                | join(per_sample_stats.map{items -> [items[0], items[1]]}, by: 0)
                | join(reference_channel, by: 0)
            )
            // Join the results summaries to the preprocessing stats
            aggregated_results = per_sample_stats
                | join(results.out.for_report, by: 0)
        } else {
            aggregated_results = per_sample_stats
        }
        
        // Collect files in a working dir, so they can be used in the final report generation 
        results_for_report = aggregated_results
            | map { items ->
                def meta = items[0]
                def stats_files = items[1..-1]
                [meta, stats_files, meta.alias]
            }
            | collectFilesInDir
            | map { _meta, dirname -> dirname }
        // if we did mapping, make the alignment sections of the report
        boolean make_mapping_report = params.mapping

        makeReport(
            metadata,
            software_versions,
            workflow_params,
            results_for_report | collect,
            make_mapping_report
        )

}

// entrypoint workflow
workflow {
    WorkflowMain.initialise(workflow, params, log)

    Pinguscript.ping_start(nextflow, workflow, params)

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
    // Check provided reference exists
    params.reference = params.reference ?: "${projectDir}/data/HG002qpMP_reference.fasta.gz"
    assert file(params.reference).exists() : "ERROR: Provided reference does not exist: ${params.reference}"

    pipeline(samples)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}