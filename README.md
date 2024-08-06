# Telomere sequencing workflow

A workflow to analyse telomere-enriched data generated using Oxford Nanopore’s Telo-Seq method.



## Introduction

`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm. The experimental protocol and sequencing are described in separate Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents. Telo-seq libraries are sequenced on Oxford Nanopore’s sequencing devices.

`wf-teloseq` has three alternative pathways to analyse Telo-Seq data: 

•	Pathway 1: Telomeric read counts and overall telomere length only.

•	Pathway 2: Chromosome arm assigned telomere read counts and telomere lengths using a matching reference. 

•	Pathway 3: Chromosome arm assigned telomere read counts and telomere lengths using a reference created de novo from the data as specific to input data.

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
[Docker](https://www.docker.com/products/docker-desktop
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
wget #TODO:
tar -xzvf wf-teloseq-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-teloseq \
	// TODO: this should use demo data--fastq './test_data/HG002_small_test.fastq.gz' \
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
| denovo | boolean | create de novo reference | If selected then the de novo guided reference is constructed by first mapping to the default or provided reference then extracting subsets of reads for clustering, consensus and polishing. This creates a reference based upon the data to map back to and separate out the telomere reads by chromosome arm | False |
| curation | boolean | add manual contigs to reference | If selected with --curation option then these user input contigs are incorporated into the reference via mapping and error correction | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| mapq | integer | mapping quality filter parameter | Mapping quality used to filter the bam file | 4 |
| read_quality | integer | read quality filter parameter | Read quality used to filter raw reads | 9 |
| enzyme_cut | string | enzyme cut site | Restriction enzyme cut site used for filtering reads that are not close to this site for the strict setting | GATATC |
| denovoRef | string | de novo assembled reference from denovo route | If provided with curatedContigs then will incorporate both sets of contigs into one reference |  |
| cov_4cluster | integer | minimum read number in clustering to produce contig, recommend 30 for 220x and 20k telomere reads and 8 for 3k telomere reads. | The minimum read number to be used for filtering after using the clustering algorithm | 7 |
| mincoverage | integer | minimum read number for coverage of reference contig, default is 20% of telomere read average for 92 chr arms | The minimum telomere coverage of chromosome arms to be taken to the final results and plots is calculated as 20% of the average chr arm coverage but can be overridden by giving a value here | -1 |
| curatedContigs | string | manually selected contigs to add to reference | If provided with denovoRef then will incorporate both sets of contigs into one reference |  |






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

`wf-teloseq` requires basecalled nanopore reads in fastq or bam format as input. If using Pathway 1, a reference genome is also required as input but not used. HG002 reference (HG002qpMP_reference.fasta.gz) in the test_data folder can be provided by the user for Pathway 1 and 3. 

We recommend basecalling you nanopore reads with the super accuracy (SUP) model during sequencing. High accuracy (HAC) basecalling will also work but will result in slightly reduced number of telomere reads. See our Know-How document for further information on this. For instructions on how to basecall your nanopore raw data, please see [dorado](https://github.com/nanoporetech/dorado). 

The recommended telomere read number input as minimum for best results would be 2000 for pathway 2 and 3000 for pathway 3.

`wf-teloseq` uses a “telomere reference” built from a full reference genome provided by the user (Pathway 2) or from the Telo-Seq data itself via a de novo reference guided approach (Pathway 3). Using a provided reference without the --denovo option will in silico lead the pipeline to find the first EcoRV enzyme cut site from the telomeric end in the sub-telomeric regions for each chromosome arm and add another 300bp. Both options of reference use will result in the chromosome arm ends oriented 5'->3' such that they all contain the telomere repeat TAACCC in the same orientation.



## Output 

`wf-teloseq` outputs telomere reads in bam format and telomere length statistics in csv format. It also produces an html report to summarise run results. If using Pathway 1, an overall sample telomere length estimate is given. Unless Pathway 2 is deactivated, the pipeline also reports a weighted average based on mean telomere lengths of chromosome arm assigned reads. This weighted average takes into account the variability of each chromosome arms read number and therefore may have less variance when comparing data sets. Only chromosome arms (or contigs) with at least 10x coverage will be reported.

Pathway 1: With the --skipmapping parameter
```
output
├── execution
│   ├── report.html
│   ├── timeline.html
│   └── trace.txt
├── reference
│   └── reference.fasta
├── Sample
│       ├── plots
│       │     └── Sample_raw_Boxplot_of_Telomere_length.pdf
│       ├── reads
│       │     └── Telomere_reads.fastq
│       └── results
│             └── Sample_raw_Per_Read_telomere_length.csv
├── params.json
├── versions.txt
└── wf-teloseq-report.html
```


Pathway 2/3: Without the --skipmapping parameter. If using --denovo then an additional de novo reference will be reported in the reference folder. This can be followed up with --curation if the user has additional curated contigs to add to the de novo reference (see https://community.nanoporetech.com/knowledge/know-how/TELO-seq on how to do this).
```
output
├── execution
│   ├── report.html
│   ├── timeline.html 
│   └── trace.txt
├── reference
│   └── reference.fasta
├── Sample
│       ├── plots   
│       │   ├── nofiltered_chr_arm_Boxplot_of_Telomere_length.pdf
│       │   ├── nofiltered_Boxplot_of_Telomere_length.pdf
│       │   ├── lowfiltered_chr_arm_Boxplot_of_Telomere_length.pdf
│       │   ├── lowfiltered_Boxplot_of_Telomere_length.pdf
│       │   ├── highfiltered_chr_arm_Boxplot_of_Telomere_length.pdf
│       │   ├── highfiltered_Boxplot_of_Telomere_length.pdf  
│       │   └── Sample_raw_Boxplot_of_Telomere_length.pdf
│       ├── reads
│       │   └── Telomere_reads.fastq
│       ├── reference
│       │   └── denovo_reference.fasta (--denovo)
│       │   └── Manual_denovo_reference.fasta (--curation)
│       ├── results        
│       │   ├── Sample_raw_Per_Read_telomere_length.csv
│       │   ├── nofiltered_Per_Read_telomere_length.csv
│       │   ├── nofiltered_chr_arm_Coverage.csv
│       │   ├── lowfiltered_Per_Read_telomere_length.csv
│       │   ├── lowfiltered_chr_arm_Coverage.csv
│       │   ├── highfiltered_Per_Read_telomere_length.csv 
│       │   └── highfiltered_chr_arm_Coverage.csv
│       └── alignments
│           ├── telomere.q4.bam
│           ├── telomere.q4.bam.bai
│           ├── highfiltered.bam
│           ├── highfiltered.bam.bai 
│           ├── lowfiltered.bam  
│           ├── lowfiltered.bam.bai
│           └── mapping_reference.fasta
├── params.json
├── versions.txt
└── wf-teloseq-report.html
```


## Different filters applied in the pipeline

| Filter | Description |
|-|-|
| Applied to all | mapq score default is 4 but user can lower or raise if needed, recommended set 0 for denovo pathway. |   
| none | no additional filters are appied so no additional reads removed. |   
| low stringency | keep only reads in which the end mapping position is 80 bp beyond last telomere motif. This is to remove short telomere only reads that would not be chromosome arm specific and also could be truncated. |
| high stringency | keep only reads in which the start mapping position is before last telomere motif identification and end mapping position is within 25 bp of cutsite with exception of cutsites beyond 45k as will get very few of these reads. This is to ensure reads span subtelomere and to limit mismapping and fragmented reads. |



## Test data

A small test dataset of reads is provided with `wf-teloseq` in "test_data" to help users test the workflow using Pathway 1 and 2. Note that this small test dataset is not big enough to test Pathway 3. It consists of a reference file and HG002 sample reads to test the installation.


## Reference genome

The telomere reference provided with the pipeline is based on the HG002 telomere-to-telomere reference genome. The pipeline requires all chromosome arms to be cut and orientated to 5'->3'. If you provide a different reference to the one in test_data folder, then the pipeline script will automatically cut the reference at enzyme cut sites +300bp where arms have telomere and extract to a new reference. Consequently, if using bam files output from this pipeline (e.g. for IGV visualisation), always use the reference output by the pipeline. 

In the supplied HG002 reference genome, the sub-telomeric sequence up to the EcoRV cut site in chromosome 13 (paternal) P arm is identical to chromosome 22 (paternal) P arm, so we only have one representative sequence in the reference provided. Although the telomere lengths cluster distinctly for these identical sequence identity arms, if two identical contigs are present there will be random mapping and a distorted telomere length estimate. The pipeline does not yet separate out these two arms based upon telomere length from the single contig provided but is planned in a future release.

Human cell lines, and individuals, have genetic variation in their sub-telomeric regions. These may impact chromosome arm assignment if differences compared to the reference used in the analysis are substantial, or make the sub-telomeric regions indistinguishable. For most reliable chromosome arm telomere length estimation, it is recommended to use a genome reference for that specific sample (e.g. cell line or individual). Overall telomere length estimation (pathway 1) is not affected by this, as it is entirely reference-free.  

`wf-teloseq` contains an experimental feature that aims to create a reference that better represents the sample. This can be called using the parameter --denovo. This approach maps telomere reads against an initial reference (HG002 by default) and then clusters reads into contigs based on differences between these reads. This feature aims to provide an improved separation of reads to chromosome arms and thus more accurate telomere length estimation.

A follow up to reference curation can be employed to change the denovo produced reference then additional contigs identified in IGV by extracting reads that represent sequences not in reference and adding these via option --curation and --curatedContigs <fasta>.



## Run time

Running a typical Telo-Seq dataset (807K) through `wf-teloseq` with human chromosome arm mapping assignment (pathway 2) takes approximately 15 minutes, and chromosome arm de novo (pathway 3) using 3000 telomere read takes 60 minutes using 16 threads and 8GB RAM. Without mapping (pathway 2), it takes less than 5 minutes. Information on typicaly Telo-Seq output can be found in our separate Know-How document but an example output report is provided in the "example_output" folder. 



## Running in Epi2me labs via Windows on a laptop 

Install the EPI2ME Desktop application [epi2me](https://labs.epi2me.io/downloads/). In the Workflows tab, select import and enter the url for this github repository (https://github.com/nanoporetech/wf-teloseq) to download the workflow. If restrictions are still in place or you have difficulty, then download the code and put in the local "workflows" folder for your epi2me installation then the workflow will be installed and available.

Pipeline outputs can be found in the publishDir location, shown in the "Parameters" box within your run analysis page.


## Example Linux commands

The pipeline is designed to be run after Dorado basecalling with input reads in (unaligned) bam or fastq files. If the files are in a folder structure i.e. 1.fastq in folder "1" and then 2.fastq in "2" folder at the same level then each sample will be run independently but presented as different tabs in the final report. `wf-teloseq` can be run with docker, and singularity.

Mapping and non-mapping analysis 

1. Pathway 1: single sample analysis without chromosome arm analysis - ~5 min
 
```nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --fastq ./test_data/ont1_test_dataset.fastq.gz --skipmapping -profile singularity ```

2. Pathway 2: single sample chromosome arm and sample analysis - ~60-180 min (8 cpu / 16 threads - 4 cpu / 8 threads)
 
```nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --fastq ./test_data/ont1_test_dataset.fastq.gz -profile singularity ```

3. Pathway 2: multi-sample analysis example
   
```nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --fastq /folder/ -profile singularity ```

4. Pathway 3: de novo reference example. If your sample contains >20K telomere reads, then setting --cov_4cluster to 30 instead of the default 8 is recommended.
   
```nextflow run main.nf --reference ./test_data/HG002qpMP_reference.fasta.gz --fastq ./test_data/ont1_test_dataset.fastq.gz --denovo --mapq 0 -profile singularity ```

5. Pathway 3b: Following from Pathway 3 with additional, user curated contigs. Uses Telomere_reads.fastq and denovo_reference5.fasta from example 4 above.
   
```nextflow main.nf --fastq Telomere_reads.fastq --curatedContigs Extra_U2_OS.fasta --denovoRef denovo_reference5.fasta --reference ./test_data/HG002qpMP_reference.fasta.gz --curation -profile singularity```


<br><br>

| Parameter | Description | options (Default if one option shown) |
|-|-|-|  
| --reference | Path to reference genome fasta file (required, even with --skipmapping) | fasta,fasta.gz |
| --fastq | Path to input reads (required) | fastq,fastq.gz |
| --bam | Path to input reads (required) | bam |
| --mapq | Minimum MAPQ quality score | 4 |
| --read_quality | Minimum read quality score | 9 |
| --cov_4cluster | Paramater for minimum number of reads for clustering if denovo or curation pathways | 8 |
| --sample_sheet | Path to sample sheet for multiple samples meta data | none |
| --skipmapping | Flag to skip mapping | False |
| --curatedContigs | Path to extracted raw contigs for addition to reference before mapping | none |
| --denovoRef | Path to denovo reference from previous run | none |
| --denovo | Flag to create reference | False |
| --curation | Flag to add to reference before mapping pipeline | False |
| --doublestranded | Flag to search both strands, note GGGTTA reads are reverse complemented to match reference orientation | False |
| --enzyme_cut | cut site of enzyme used, 1 only | GATATC |
| -profile | Environment profile (required) | docker,singularity |



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



