#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { getParams } from './lib/common'

include { align_and_process } from './modules/local/alignment'
include { collectFilesInDir } from './modules/local/common'
include { process_reads } from './modules/local/preprocess'

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
    python -c "import edlib; import importlib.metadata; print(f'edlib,{importlib.metadata.version(\\"edlib\\")}')" >> versions.txt
    python -c "import scipy; print(f'scipy,{scipy.__version__}')" >> versions.txt
    samtools --version | grep samtools | head -1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | awk '{print "minimap2," \$1}' >> versions.txt
    """
}

// Make report html file using staged file dir and other files
process makeReport {
    label "wf_common"
    cpus   1
    memory '2 GB'
    publishDir (
        params.out_dir,
        mode: "move",
    )
    input:
        val meta_array
        path "versions/*"
        path "params.json"
        path "data/*"
    output:
        path "wf-teloseq-*.html"
    script:
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
            --workflow-version v${workflow.manifest.version}
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
            .filter { _meta, sub_telomeric_fastq_file, _length_sum, _stats_dir ->
                sub_telomeric_fastq_file.size() > 0
            }
            .ifEmpty {
                throw new Exception("No valid samples found. Exiting workflow.")
            }

    emit:
        metadata = valid_samples_with_stats.map { items -> items[0] }.collect()
        valid_samples_with_stats
}


workflow pipeline {
    take:
        samples

    main:
        // Get the workflow_params and software versions for report.
        workflow_params = getParams()
        software_versions = getVersions()


        preprocessing = generic_preprocessing(samples)
        
        // Extract outputs from preprocessing
        metadata = preprocessing.metadata
        per_sample_stats = preprocessing.valid_samples_with_stats

        // Define default reference, replaced by user option if supplied.
        def default_ref = params.reference ? file(params.reference, checkIfExists: true) : file("$projectDir/data/HG002qpMP_reference.fasta.gz", checkIfExists: true)

        if (params.mapping) {
            //map filtered telomere reads to genome and filter using mapq (default=4)
            aligned_outputs = align_and_process(per_sample_stats.map{items -> 
                def meta = items[0]
                def reads = items[1]
                // Collect statistics files from fastcat and preprocessing
                def stats_files = items[2..-1]
                // Use default reference if user hasn't set one 
                if (!meta.containsKey('reference')) {
                    meta.reference = default_ref
                }
                [meta, reads, file(meta.reference), stats_files]
            })
            analysis_files_for_collection = aligned_outputs.alignment_stats
        } else {
            analysis_files_for_collection = per_sample_stats
        }
        // Collect files in a working dir, so they can be used in the final report generation 
        results_for_report = analysis_files_for_collection
            | map { items ->
                def meta = items[0]
                // The fastcat and preprocessing stats files are nested, so flatten them out
                def stats_files = items[1..-1].flatten()
                [meta, stats_files, meta.alias]
            }
            | collectFilesInDir
            | map { _meta, dirname -> dirname }

        makeReport(
            metadata,
            software_versions,
            workflow_params,
            results_for_report | collect,
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