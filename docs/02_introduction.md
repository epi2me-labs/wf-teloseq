`wf-teloseq` is a Nextflow pipeline to provide telomere length estimates from data generated with the Telo-Seq protocol. 

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm.
`wf-teloseq` currently defaults to an alignment based analysis.
By default the workflow will use a suitable reference which has been created from the telomeric regions of the phased telomere-to-telomere (T2T) HG002 human reference genome, but users may optionally provide their own reference genome.
Please note that a T2T reference genome matching the sample that has been sequenced is required for correct assignment of the reads to chromosome arms.
If no reference genome is available, a reference-less analysis can be performed, which simply provides a global telomere length estimate based on the sequenced reads without chromosome arm assignment.

**NOTE**
In order to run the analysis, data must first be base-called and demultiplexed with Dorado using Telo-Seq specific configuration files.
A bash script to perform this base-calling is available in the `wf-teloseq` repository under `bin/basecalling.sh`.
