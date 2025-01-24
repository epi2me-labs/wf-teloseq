# Telomere sequencing workflow

A workflow to analyse telomere-enriched data generated using Oxford Nanopore’s Telo-Seq method.



## Introduction

`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm. The experimental protocol and sequencing are described in separate Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents. Telo-seq libraries are sequenced on Oxford Nanopore’s sequencing devices.

`wf-teloseq` currently supports two alternative pathways to analyse Telo-Seq data: 

•	Pathway 1: Overall telomeric read counts and telomere length only (i.e. this combines all reads, regardless of which chromosome they originated from).

•	Pathway 2: Using a matching reference, determine telomere read length and count for each chromosome arm individually.

## Note

`wf-teloseq` is in an early development phase and may undergo changes, improvements, and feature additions in the future. 



## Compute requirements

Recommended requirements:

+ CPUs = 16
+ Memory = 32GB

Minimum requirements:

+ CPUs = 8
+ Memory = 16GB

Approximate run time: 5/15/30 minutes per sample for pathway 1/2/3

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

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

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
| doublestranded | boolean | double stranded protocol | If selected then the telomere reads are identified on both strands, then the second strand is reverse complemented so that all reads are orientated telomere first and pipeline continues as normal but with the other strand reads used. | False |
| reference | string | Reference genome of the sequenced sample. | Reference genome of the sequenced material in fasta format. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### TeloSeq Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| skipmapping | boolean | Skip mapping step for just sample only telomere length | If selected then the workflow will not run the mapping step but measure telomere length just on the unmapped telomere identified reads. | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| adapters_to_trim | string | Comma-delimited list of sequencing adapters to be trimmed |  | CACCCTAACCCTAACCCTAACC,CAACCCTAACCCTAACCCTAAC,CCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCT,CTAACCCTAACCCTAACCCTAA,CCCCTAACCCTAACCCTAACCC |
| filter_error_motifs | string | Comma-delimited list of telomere error motifs for filtering | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in a telomere read, it is dropped from the analysis. | GTATAG,CGCGCGCG,CCACCG,AGCGACAG,ATAAGT,CCTCGTCC,TATAGT,AGTACT,GAGTCC,TATAGT,TATACA,TGGTCC,CTCTCCTCT |
| filter_error_motifs_max_count | integer | Number of erroneous telomere k-mers to remove read. | If more than this number of telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in the telomeric region of a read, it is dropped from the analysis. | 5 |
| filter_error_motifs_window_size | integer | Size of window for filtering based on telomere error motifs. | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of this length in the telomeric region of a read, it is dropped from the analysis. | 500 |
| mapq | integer | mapping quality filter parameter | Mapping quality used to filter the bam file | 4 |
| read_quality | integer | read quality filter parameter | Read quality used to filter raw reads | 9 |
| restriction_site | string | enzyme cut site | Restriction enzyme cut site used for filtering reads that are not close to this site for the strict setting | GATATC |
| cov_4cluster | integer | minimum read number in clustering to produce contig, recommend 30 for 220x and 20k telomere reads and 8 for 3k telomere reads. | The minimum read number to be used for filtering after using the clustering algorithm | 7 |
| mincoverage | integer | minimum read number for coverage of reference contig, default is 20% of telomere read average for 92 chr arms | The minimum telomere coverage of chromosome arms to be taken to the final results and plots is calculated as 20% of the average chr arm coverage but can be overridden by giving a value here | -1 |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-teloseq-report.html | Report for all samples | aggregated |
| Per file read stats | ./fastq_ingress_results/reads/fastcat_stats/per-file-stats.tsv | A TSV with per file read stats, including all samples. | aggregated |
| Per read stats | ./fastq_ingress_results/reads/fastcat_stats/per-read-stats.tsv | A TSV with per read stats, including all samples. | aggregated |
| Run ID's | ./fastq_ingress_results/reads/fastcat_stats/run_ids | List of run ID's present in reads. | aggregated |
| Meta map json | ./fastq_ingress_results/reads/metamap.json | Meta data used in workflow presented in a JSON. | aggregated |
| Concatenated sequence data | ./fastq_ingress_results/reads/{{ alias }}.fastq.gz | Per sample reads concatenated in to one fastq file. | per-sample |




## Pipeline overview

## Input

`wf-teloseq` requires basecalled nanopore reads in fastq or bam format as input.
For best results, we recommend to use at least 2000 telomere reads per sample.
Unless per-chromosome analysis is disabled, a reference is also required.
This can be a telomere-to-telomere assembly or a FASTA file containing telomeres only.
Per default, a telomere reference constructed from HG002 (`data/HG002qpMP_reference.fasta.gz`) is used.

We recommend basecalling with super accuracy (SUP) models.
High accuracy (HAC) basecalling will also work but will result in slightly reduced number of telomere reads.
See our [Know-How document](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) for further information on this.

The workflow will search both ends of each reference sequence for telomeric repeats (of `TAACCC`) and the first occurrence of the restriction enzyme motif provided by `--restriction_site` (per default `GATATC` for EcoRV).
It will then
* extract the telomere until `--ref_add_margin` after the restriction site
* attempt to reduce bias in the later alignment step caused by length differences and variants in the telomeric sequences by
    * trimming the telomere sequence (keeping only the sub-telomere)
    * extending the sequence with a "synthetic" telomere sequence of 4000 repeats of `TAACCC`

The reads are then aligned against the reference + filtered (see below for details) before determining telomere length + counts.


## Output 

`wf-teloseq` outputs aligned telomere reads in bam format and telomere length statistics in csv format. It also produces an html report to summarise run results. Unless reference analysis is disabled, an overall telomere length estimate is given for each sample. Otherwise, the pipeline also reports the mean telomere lengths of reads assigned to individual chromosome arms. Only chromosome arms (or contigs) with at least 10x coverage will be reported.


## Different filters applied in the pipeline

| Filter | Description |
|-|-|
| Applied to all | mapq score default is 4 but user can lower or raise if needed |   
| none | no additional filters are appied so no additional reads removed. |   
| low stringency | keep only reads in which the end mapping position is 80 bp beyond last telomere motif. This is to remove short telomere only reads that would not be chromosome arm specific and also could be truncated. |
| high stringency | keep only reads in which the start mapping position is before last telomere motif identification and end mapping position is within 25 bp of cutsite with exception of cutsites beyond 45k as will get very few of these reads. This is to ensure reads span subtelomere and to limit mismapping and fragmented reads. |



## Test data

A small test dataset of reads is provided with `wf-teloseq` in "test_data" to help users test the workflow using Pathway 1 and 2. Note that this small test dataset is not big enough to test Pathway 3. It consists of a reference file and HG002 sample reads to test the installation.


## Reference genome

The telomere reference provided with the pipeline is based on the HG002 telomere-to-telomere reference genome. The pipeline requires all chromosome arms to be cut and orientated to 5'->3'. If you provide a different reference to the one in test_data folder, then the pipeline script will automatically cut the reference at enzyme cut sites +300bp where arms have telomere and extract to a new reference. Consequently, if using bam files output from this pipeline (e.g. for IGV visualisation), always use the reference output by the pipeline. 

In the supplied HG002 reference genome, the sub-telomeric sequence up to the EcoRV cut site in chromosome 13 (paternal) P arm is identical to chromosome 22 (paternal) P arm, so we only have one representative sequence in the reference provided. Although the telomere lengths cluster distinctly for these identical sequence identity arms, if two identical contigs are present there will be random mapping and a distorted telomere length estimate. The pipeline does not yet separate out these two arms based upon telomere length from the single contig provided but is planned in a future release.

Human cell lines, and individuals, have genetic variation in their sub-telomeric regions. These may impact chromosome arm assignment if differences compared to the reference used in the analysis are substantial, or make the sub-telomeric regions indistinguishable. For most reliable chromosome arm telomere length estimation, it is recommended to use a genome reference for that specific sample (e.g. cell line or individual). Overall telomere length estimation (pathway 1) is not affected by this, as it is entirely reference-free.  


## Run time

Running a typical Telo-Seq dataset (807K) through `wf-teloseq` with human chromosome arm mapping assignment takes approximately 15 minutes. When skipping the mapping stage (`--skipmappnig`), it takes less than 5 minutes.


## Running in Epi2me labs via Windows on a laptop 

Install the EPI2ME Desktop application [epi2me](https://labs.epi2me.io/downloads/). In the Workflows tab, select import and enter the url for this github repository (https://github.com/nanoporetech/wf-teloseq) to download the workflow. If restrictions are still in place or you have difficulty, then download the code and put in the local "workflows" folder for your epi2me installation then the workflow will be installed and available.

Pipeline outputs can be found in the publishDir location, shown in the "Parameters" box within your run analysis page.




## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

## Additional

Telo-Seq has only been tested on human data, but we expect it to work well on species with similar telomere repeat sequences.

Telomere sequences are highly repetitive. Sometimes basecalling errors cause dense mis-basecalled repeats to occur in Telo-seq reads. These get soft-clipped off by the aligner in downstream analysis, which could result in underestimation of telomere lengths. `wf-teloseq` identifies such reads and excludes them from telomere length estimation.


# Acknowledgements

This project uses code in the "telomerewindowV1.py" script from the following source:

- **Original Author:** Ramin Kahidi
- **Original Repository:** https://github.com/GreiderLab/TeloBP

The code used is licensed under the MIT License, which can be found in the original repository and the header of the script.

Publication: Karimian K, Groot A, Huso V, Kahidi R, Tan KT, Sholes S, Keener R, McDyer JF, Alder JK, Li H, Rechtsteiner A, Greider CW. Human telomere length is chromosome specific and conserved across individuals. 2024 Jan 13:2023.12.21.572870. doi: 10.1101/2023.12.21.572870. Update in: Science. 2024 May 3;384(6695):533-539. doi: 10.1126/science.ado0431. PMID: 38187739; PMCID: PMC10769321.


# Epi2Me

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



