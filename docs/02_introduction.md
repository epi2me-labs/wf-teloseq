`wf-teloseq` is a bioinformatics pipeline to analyse nanopore telomere sequencing (Telo-Seq) data. It is implemented in Nextflow.

Telo-Seq aims to measure telomere length accurately and assign each telomere to a chromosome arm. The experimental protocol and sequencing are described in separate Telo-Seq protocol and [Know-How](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) documents. Telo-seq libraries are sequenced on Oxford Nanopore’s sequencing devices.

`wf-teloseq` has three alternative pathways to analyse Telo-Seq data: 

•	Pathway 1: Telomeric read counts and overall telomere length only.

•	Pathway 2: Chromosome arm assigned telomere read counts and telomere lengths using a matching reference. 

•	Pathway 3: Chromosome arm assigned telomere read counts and telomere lengths using a reference created de novo from the data as specific to input data.

## Note

`wf-teloseq` is in an early development phase and may undergo changes, improvements, and feature additions in the future. 