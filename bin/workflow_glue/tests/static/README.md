### Test data for workflow glue tests

#### test_telomere_boundary.fastq.gz
Iterated the wf-teloseq test data set, using `process_reads.py`.
Saved the first three reads which hit a certain BoundaryFinder variant into this FASTQ.

#### test_main_alignment.bam
Reads from the test-dataset with one added from the release dataset barcode 01, which QC's as `TooCloseStart`.

#### short_subtelomere.bam
Reads which align to Chr5 Paternal P. This contig arm has a short subtelomere section, which can be misidentified as `TooCloseEnd`.

#### test_main.bam
Took 10 reads from Barcode01 of the release dataset, from every QC mode assigned.
Merged these reads with the contents of `test_telomere_boundary.bam`.
Hits every QC mode, used to test `process_reads`.