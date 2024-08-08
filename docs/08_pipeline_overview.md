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
