### Test data for workflow glue tests

#### test_ref_reads.bam
Calculated the boundary of telomere repeats for each contig in the reference contained in wf-teloseq/data/HG002qpMP_reference.fasta.gz.
Created reads which are boundary -1000, boundary + 2000.
Generated between 0 and 10 reads for each contig, mapped these back to wf-teloseq/data/HG002qpMP_reference.fasta.gz.
Added in the TL tag for these reads as calculated by `find_telomere_boundary` in `process_reads.py`.
Used for tests of `process_alignments.py`.

#### test_telomere_boundary.fastq.gz
Iterated the wf-teloseq test data set, using `process_reads.py`.
Saved the first three reads which hit a certain BoundaryFinder variant into this FASTQ.