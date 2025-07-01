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
        mode: "move",
    )
    input:
        val meta_array
        path "versions/*"
        path "versions/*"
        path "params.json"
        path unaligned_stats, stageAs: "unaligned_stats/*"
        path kde_stats, stageAs: "kde_stats/*"
        path fastcat_stats_dirs, stageAs: "fastcat_stats/staged_*"
        path qc_stats, stageAs: "qc_stats/*"
        path contig_stats, stageAs: "contig_stats/*"
        path boxplot_stats, stageAs: "boxplot_stats/*"
    output:
        path "${report_name}"
    script:
        report_name = "wf-teloseq-report.html"  // nodef: used by output
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
            --unaligned-stats-dir unaligned_stats \\
            --kde-stats-dir kde_stats \\
            --fastcat-stats-dir fastcat_stats \\
            --qc-stats-dir qc_stats \\
            --contig-stats-dir contig_stats \\
            --boxplot-stats-dir boxplot_stats
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
            .map { meta, sub_telomeric_fastq_file, length_sum, kde_sum, stats_dir ->
                if (sub_telomeric_fastq_file.size() == 0) {
                     log.warn "SAMPLE: ${meta.alias} contains no valid sequences, removed from analysis..."
                }
                else {
                    [meta, sub_telomeric_fastq_file, length_sum, kde_sum, stats_dir]
                }
            }
            .filter { _meta, sub_telomeric_fastq_file, _length_sum, _kde_sum, _stats_dir ->
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
            // map filtered telomere reads to genome
            aligned_outputs = align_and_process(per_sample_stats.map{items -> 
                def meta = items[0]
                def reads = items[1]
                // Collect statistics files from fastcat and preprocessing
                def stats_files = items[2..-1]
                // Use default reference if user hasn't set one 
                if (!meta.containsKey('reference')) {
                    meta.reference = params.reference
                }
                [meta, reads, file(meta.reference), stats_files]
            })
            analysis_files_for_collection = aligned_outputs.alignment_stats
        } else {
            // map to get only stats files out -> matches to alignment output
            analysis_files_for_collection = per_sample_stats.map{ meta, _sub_telomeric_fastq_file, length_sum, kde_stats, stats_dir -> 
                [meta, length_sum, kde_stats, stats_dir]
            }
        }
        // Sort files from output tuple into individual channels so they can be collected and staged in report working directory.
        // We use the whole `items` tuple as we don't know how many files are contained.
        // If alignment is performed it will contain extra stats files.
        // Finally, we use the order of samples in the collected metadata_ch to match up to the order of the fastcat_stats directories staged in `makeReport`.
        (
            metadata_ch,
            unaligned_ch,
            kde_ch,
            fastcat_stats_ch,
            boxplot_stats_ch,
            qc_stats_ch,
            contig_stats_ch
        ) = analysis_files_for_collection.multiMap { items ->
            def meta = items[0]

            // The fastcat and preprocessing stats files are nested, so flatten everything out
            def stats_files = items[1..-1].flatten()
            // These files are only output if alignment is performed, set a default.
            def boxplot_stats = OPTIONAL_FILE
            def qc_stats = OPTIONAL_FILE
            def contig_summary_stats = OPTIONAL_FILE
            //  Alignment was performed, set stats files paths for report.
            if (!params.skip_mapping) {
                  boxplot_stats = stats_files[4]
                  qc_stats = stats_files[5]
                  contig_summary_stats = stats_files[6]
            }

            // Return file elements
            metadata_ch: meta
            unaligned_stats_ch: stats_files[0]
            kde_ch: stats_files[1]
            fastcat_stats_ch: stats_files[2]
            boxplot_stats_ch: boxplot_stats
            qc_stats_ch: qc_stats
            contig_summary_stats_ch: contig_summary_stats
        }

        makeReport(
            metadata_ch.toList(),
            software_versions,
            software_version_common,
            workflow_params,
            unaligned_ch | collect,
            kde_ch | collect,
            fastcat_stats_ch.toList(),
            qc_stats_ch | collect,
            contig_stats_ch | collect,
            boxplot_stats_ch | collect
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
