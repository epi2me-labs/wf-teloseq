# Telomere sequencing workflow

A workflow to analyse telomere-enriched data generated using Oxford Nanopore’s Telo-Seq method.



## Introduction

`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm.
The experimental protocol and sequencing are described in the Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents.
Telo-seq libraries are sequenced on Oxford Nanopore Technologies' sequencing devices.

`wf-teloseq` currently defaults to an alignment based analysis.
Users can provide a reference, or use a default which has been created from the telomeric regions of the phased telomere-to-telomere HG002 human reference genome.
Alternatively, a reference-less "bulk" analysis can be performed, which simply provides a combined length estimation for all telomeres in each sample, for use in cases where a suitable reference is not available. 

> In order to run the analysis data must have been basecalled and demultiplexed using Telo-Seq specific basecalling.
> A bash script to perform this basecalling is available in the wf-teloseq repository under `bin/basecalling.sh`. 




## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 32GB

Minimum requirements:

+ CPUs = 8
+ Memory = 16GB

Approximate run time: Approx 5 minutes per sample, or about 1 minute if not performing alignment.

ARM processor support: False




## Install and run


These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop)
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-teloseq --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-teloseq
```

A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-teloseq/wf-teloseq-demo.tar.gz
tar -xzvf wf-teloseq-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-teloseq \
	--bam 'wf-teloseq-demo/teloseq_example' \
	--reference 'wf-teloseq-demo/HG002qpMP_reference.fasta.gz' \
	-profile standard
```

For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/




## Related protocols

<!---Hyperlinks to any related protocols that are directly related to this workflow, check the community for any such protocols.--->

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices.

Find related protocols in the [Nanopore community](https://community.nanoporetech.com/docs/).



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. The sample sheet can include a "reference" column to assign each sample a specific reference to map to.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| reference | string | Reference genome of the sequenced sample, if not specified a human derived telomere reference set will be used. |  |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .FASTQ(.gz) file or directory of .FASTQ(.gz) files. |  |  |


### TeloSeq Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| skip_mapping | boolean | Perform alignment to assign haplotypes to telomeric reads. | Use `--skip_mapping` if there is no suitable reference available. Only a bulk estimate of telomere lengths per sample will be calculated. | False |
| alignment_threads | integer | Set max number of threads to use for alignment. |  | 6 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| min_length | integer | Minimum read length for filtering. | Used in the initial filtering of reads into the workflow. Reads with a sequence length less than this will be removed prior to analysis, and will not be present in the output files. This removes all reads which are either noise, artifacts or partial, and which cannot be full telomeric sequences. Setting this number too high runs the risk of removing short telomeres. | 100 |
| read_quality | integer | Minimum read Q score for filtering. | Reads with a mean quality score lower than this value will be filtered out prior to analysis, and will not be present in the output files. | 9 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-teloseq-report.html | Report for all samples. | aggregated |
| Tool versions | versions.txt | A CSV with per row tool and version. | aggregated |
| Parameters from workflow | params.json | A json of all parameters selected in workflow. | aggregated |
| Unaligned, filtered and tagged sequences. | {{alias}}/unaligned_data/{{alias}}_filtered_telomeric.fastq | These sequences have been tagged (valid SAM format tags) on whether or not they passed filtering (qc:Z), and if detected, also tagged with the Telomere repeat boundary coordinates (tl:I). | per-sample |
| Summary metrics about detected telomere lengths within the sample. | {{alias}}/unaligned_data/{{alias}}_telomere_unaligned_metrics.tsv | Aggregated summary metrics about reads which have passed all filtering, with a detected telomere repeat boundary. | per-sample |
| Aligned, filtered and tagged sequences. | {{alias}}/aligned_data/alignment_filtered_teloseqs.bam | Contains the sequences from the processed_fastq aligned to the provided reference. All sequences which passed initial read length (`--min_length`) and quality (`--read_quality`) filtering will be present, including unmapped. Only produced if alignment is performed. | per-sample |
| Accompanying index for the aligned BAM. | {{alias}}/aligned_data/alignment_filtered_teloseqs.bam.csi | Coordinate-sorted index file for the aligned data BAM. | per-sample |
| Summary metrics about aligned telomere lengths within the sample. | {{alias}}/aligned_data/{{alias}}_telomere_aligned_metrics.tsv | Aggregated summary metrics about detected telomere lengths within the sample, after alignment. Only reads which have primary alignments to the reference are considered. | per-sample |
| Summary metrics about aligned telomere lengths grouped by target contig. | {{alias}}/aligned_data/{{alias}}_contig_telomere_aligned_metrics.tsv | Aggregated summary metrics about detected telomere lengths within the sample, after alignment. The reads are grouped by target contig, and only reads which have primary alignments to the reference are considered. | per-sample |
| Summary metrics about the filtering status of each read. | {{alias}}/aligned_data/{{alias}}_qc_modes_metrics.tsv | Aggregated summary metrics of the filtering status of each read. Metrics displayed for each filtering status include count, mean Q score, length, and alignment identity (where applicable). | per-sample |




## Pipeline overview

The workflow is composed of two steps. 
Firstly, reads undergo initial basic length and quality control filtering, followed by analysis to determine the telomeric boundary, and filtering to remove reads where the boundary is likely incorrect.
If alignment is enabled (default), all reads (irrespective of boundary determination) are then aligned to a reference. Reads with primary alignments with good telomere boundaries are then aggregated, before statistics about the estimated lengths of each contigs telomeresß are generated. 

## 1. Input and Sequence preparation.

`wf-teloseq` requires basecalled nanopore reads in FASTQ or BAM format as input.
We recommend basecalling with super accuracy (SUP) models.
High accuracy (HAC) basecalling will also work but will result in a slightly reduced number of telomeric reads of sufficient quality for analysis.
For best results, we recommend a minimum of 1000 telomere reads per sample for alignment based analyses.
If a reference is not provided via the `--reference` or `--sample_sheet` parameter, then a telomere reference constructed from HG002 (`data/HG002qpMP_reference.fasta.gz`) is used by default for alignment.

See our [Know-How document](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) for further information.


## 2. Read filtering and tagging
Filters are listed in the order that they are applied.
The first two filters are tunable with the params `--min_length` and `--read_quality`. 
Other filters are currently not tunable via parameters.

Reads which fail a check are tagged with a `qc:Z:<Tag>` tag, in either the output FASTQ if not performing alignment or the output BAM.
Tags for corresponding failure modes are listed below.
Reads which pass all filtering checks are tagged with `qc:Z:Good`.

| Filter               | Description                                                                 | Tag               |
|----------------------|---------------------------------------------------------------------------|------------------|
| **Minimum read length** | Reads under 100 bases long are discarded. Removed entirely.             | N/A              |
| **Minimum Q score**  | Reads with a Mean Q score <9 are discarded. Removed entirely.            | N/A              |
| **Too Short**        | The read was too short (less than 160 bases). Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooShort        |
| **Too Few Repeats**  | There were fewer than 20 telomeric repeat motifs across the entire read. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooFewRepeats   |
| **Start Not Repeats** | The first 30% of the read is not 80% repeats. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | StartNotRepeats |
| **Too Close End**    | The telomeric boundary is too close (within 80 bases) of the end of the read. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooCloseEnd     |
| **Low Sub-Telo Qual** | The mean basecall Q score of the region after the boundary is below a default value of 9. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | LowSubTeloQual  |
| **Too Errorful**     | A large number of known basecall error motifs has been observed in the subtelomere. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooErrorful     |

The following filter is only applied if alignment is performed:
| Filter               | Description                                                                 | Tag               |
|----------------------|---------------------------------------------------------------------------|------------------|
| **Bad Alignment**    | Query read has a low Gap compressed identity to the reference. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | BadAlign        |

All reads have a `tl:i` tag set, which represents the estimated telomere length. 
For reads which a boundary CANNOT be determined, the value of this tag is set to -1.

It is after this stage that a tagged FASTQ file is available per sample, with the tags present in the header.

## 3. Alignment.
By default, ALL reads are aligned by minimap2 (excepting reads which fail the initial minimum read length and Q score filters) after initial processing. 
Alignment can be skipped by setting `--mapping false` when running the pipeline.

### Default reference genome

The telomere reference provided with the pipeline is based on the [HG002 telomere-to-telomere reference genome](https://github.com/marbl/hg002).
Note that in this reference, the sub-telomeric sequence up to the EcoRV cut site in chromosome 13 (paternal) P arm is identical to chromosome 22 (paternal) P arm, so there is only one representative sequence present in the file, with a sequence name of `chr13_22_PATERNAL_P`.
The pipeline does not separate out these two arms based upon telomere length from the single contig provided, which means the estimated telomere length assigned to this contig is likely inaccurate.

### Processing of alignments
After alignment, the alignments are used to aggregate stats based on the reference contigs.
One final filtering step is first performed here, which is based on the Gap compressed identity between the query sequence and the target reference, as calculated by minimap2.
Any read with an identity of less than 0.8 is excluded from the estimated telomere length statistics.
Only primary alignments are used in this aggregation, and only those from reads which have the `qc:Z` tag set with a value of `Good`, indicating they passed all initial filtering steps.

After this stage a tagged BAM file containing both mapped and unmapped reads is produced and output per sample.




## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-teloseq/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

## Additional

Telo-Seq has only been tested on human data, but we expect it to work well on species with similar telomere repeat sequences.

Telomere sequences are highly repetitive. Sometimes basecalling errors cause dense mis-basecalled repeats to occur in Telo-seq reads. These get soft-clipped off by the aligner in downstream analysis, which could result in underestimation of telomere lengths. `wf-teloseq` identifies such reads and excludes them from telomere length estimation.

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



