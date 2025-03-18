`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm. The experimental protocol and sequencing are described in separate Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents. Telo-seq libraries are sequenced on Oxford Nanopore’s sequencing devices.

`wf-teloseq` currently defaults to an alignment based analysis. Users can provide a reference, or use a default which has been created from the telomeric regions of the T2T, haplotyped HG002 human reference genome.
Alternatively, a reference-less "bulk" analysis can be performed, which simply provides a combined length estimation for all telomeres in each sample, for use in cases where a suitable reference is not available. 


## Note

`wf-teloseq` is in an early development phase and may undergo changes, improvements, and feature additions in the future. 