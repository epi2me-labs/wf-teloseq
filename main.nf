#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { align_and_process } from './modules/local/alignment'
include { process_reads } from './modules/local/preprocess'

process getVersions {
    label "wf_teloseq"
    cpus 1
    memory '2 GB'
    output: path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas as pd; print(f'pandas,{pd.__version__}')" >> versions.txt
    python -c "import numpy as np; print(f'numpy,{np.__version__}')" >> versions.txt
    python -c "import edlib; import importlib.metadata; print(f'edlib,{importlib.metadata.version(\\"edlib\\")}')" >> versions.txt
    python -c "import ezcharts as ez; print(f'ezcharts (analysis),{ez.__version__}')" >> versions_common.txt
    python -c "import scipy; print(f'scipy,{scipy.__version__}')" >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | awk '{print "minimap2," \$1}' >> versions.txt
    """
}

process getVersionsCommon {
    label "wf_common"
    cpus 1
    memory '2 GB'
    output: path "versions_common.txt"
    script:
    """
    python -c "import ezcharts as ez; print(f'ezcharts (report),{ez.__version__}')" >> versions_common.txt
    """
}

// Make report html file using staged file dir and other files
process makeReport {
    label "wf_common"
    cpus   1
    memory '2 GB'
    publishDir (
        params.out_dir,
        mode: "copy",
        pattern: "${report_name}"
    )
    input:
        val meta_array
        path "versions/*"
        path "versions/*"
        path "params.json"
        path(stats, stageAs: "staged_*/stats")
        path(fastcat_stats, stageAs: "staged_*/fastcat_stats")
    output:
        path "${report_name}"
    script:
        report_name = "wf-teloseq-report.html"  // nodef: used by output and publishDir
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
            --report-stats-dir stats
        """
}

workflow generic_preprocessing {
    take:
        samples
    
    main:
        // Filter reads down to high quality telomere repeat spanning reads. Errors if no reads left after processing all samples
        // Saves some length based statistics about sample read datasets
        filtered_samples = process_reads(samples)

        valid_samples_with_stats = filtered_samples
            .map { meta, sub_telomeric_fastq_file, fastcat_stats_dir, stats_files_dir ->
                if (!meta.containsKey('reference')) {
                    meta.reference = params.reference
                }
                if (sub_telomeric_fastq_file.size() == 0) {
                     log.warn "SAMPLE: ${meta.alias} contains no valid sequences, removed from analysis..."
                }
                else {
                    [meta, sub_telomeric_fastq_file, fastcat_stats_dir, stats_files_dir]
                }
            }
            .filter { _meta, sub_telomeric_fastq_file, _fastcat_stats_dir, _stats_files_dir ->
                sub_telomeric_fastq_file.size() > 0
            }
            .ifEmpty {
                throw new Exception("No valid samples found. Exiting workflow.")
            }

    emit:
        valid_samples_with_stats
}


workflow pipeline {
    take:
        samples

    main:

        def OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")
        // Get the workflow_params and software versions for report.
        workflow_params = getParams()
        software_versions = getVersions()
        software_version_common = getVersionsCommon()

        preprocessing = generic_preprocessing(samples)
        
        // Extract output from preprocessing
        per_sample_stats = preprocessing.valid_samples_with_stats

        if (!params.skip_mapping) {
            per_sample_files = per_sample_stats.multiMap {
                meta, reads, fastcat_stats_dir, stats_files_dir ->
                analysis_files: [meta, reads, file(meta.reference), stats_files_dir]
                fastcat_stats_dirs: [meta, fastcat_stats_dir]

            }
            // map filtered telomere reads to genome
            aligned_outputs = align_and_process(per_sample_files.analysis_files)
            // Rejoin fastcat stats onto aligned output amd reorder tuple to match structure
            // if no alignment is performed
            analysis_files_for_collection = aligned_outputs.stats.join(
                per_sample_files.fastcat_stats_dirs, by: 0
            ).map {
                meta, stats_files_dir, fastcat_stats_dir ->
                [meta, fastcat_stats_dir, stats_files_dir]
            }
        } else {
            // map to get only stats files output -> matches same tuple structure as alignment output
            analysis_files_for_collection = per_sample_stats.map { meta, _sub_telomeric_fastq_file, fastcat_stats_dir, stats_files_dir -> 
                [meta, fastcat_stats_dir, stats_files_dir]
            }
        }
        // Split out the metadata, fastcat stats directories and stats directories into seperate channels under the 
        // report_data umbrella
        // Use the order of sample elements in the `metadata` tuple to match up to the
        // order of the `fastcat_stats` and `stats` directories after staging in `makeReport`.
        report_data = analysis_files_for_collection.multiMap { meta, fastcat_stats, stats ->
            metadata: meta
            fastcat_stats: fastcat_stats
            stats: stats
        }

        makeReport(
            report_data.metadata.toList(),
            software_versions,
            software_version_common,
            workflow_params,
            report_data.stats.toList(),
            report_data.fastcat_stats.toList()
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
    
    // use default reference if not provided, this is like this because a reference
    // isn't strictly required when skip_mapping is set, and awkwardness in validating
    // schema and params
    if (!params.reference) {
        CWUtil.mutateParam(params, "reference", "${projectDir}/data/HG002qpMP_reference.fasta.gz")
    }
    assert file(params.reference).exists() : "ERROR: Provided reference does not exist: ${params.reference}"

    pipeline(samples)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)
}
