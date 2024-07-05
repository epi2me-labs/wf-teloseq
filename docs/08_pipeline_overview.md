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