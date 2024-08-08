`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm. The experimental protocol and sequencing are described in separate Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents. Telo-seq libraries are sequenced on Oxford Nanopore’s sequencing devices.

`wf-teloseq` currently supports two alternative pathways to analyse Telo-Seq data: 

•	Pathway 1: Overall telomeric read counts and telomere length only (i.e. this combines all reads, regardless of which chromosome they originated from).

•	Pathway 2: Using a matching reference, determine telomere read length and count for each chromosome arm individually.

## Note

`wf-teloseq` is in an early development phase and may undergo changes, improvements, and feature additions in the future. 