`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm.
The experimental protocol and sequencing are described in the Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents.
Telo-seq libraries are sequenced on Oxford Nanopore Technologies' sequencing devices.

`wf-teloseq` currently defaults to an alignment based analysis.
Users can provide a reference, or use a default which has been created from the telomeric regions of the phased telomere-to-telomere HG002 human reference genome.
Alternatively, a reference-less "bulk" analysis can be performed, which simply provides a combined length estimation for all telomeres in each sample, for use in cases where a suitable reference is not available. 

> In order to run the analysis data must have been basecalled and demultiplexed using Telo-Seq specific basecalling.
> A bash script to perform this basecalling is available in the wf-teloseq repository under `bin/basecalling.sh`. 
