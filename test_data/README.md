
# TEST DATA
This data is used the CI pipeline intgration tests.

## Files

- `sample_sheet_different_ref.csv`: An example sample sheet, using the additional reference column (detailed nowhere :fire:). Provides two samples, each to be mapped to a different reference.
- `sample_sheet_same_ref.csv`: An example sample sheet, using the additional reference column (detailed nowhere :fire:). Provides two samples, each to be mapped to the same reference.
- `HG002_small_test.fastq.gz`: Small subset of 1000 filtered telomere reads to test mapping and unmapped telomere estimation pathways.
- `YAOqpMP_reference.fa.gz`: Chinese public human genome with telomere contigs extracted and trimmed to restriction enzyme cutting location, example of an alternative reference.
- `non_telomeric_1000.bam`: 1000 non telomeric human sequence reads. These reads are the first 1000 reads from the `demo.bam` file provided in the test data for `wf-human-variation`.
- `expected_out/expected_default_output.csv`: The output of the telomeric assignment step as of 026b838ec942a3cad5de1f69f99fc7e0ccd32994.
- `samples/barcode*/*.bam`: Test dataset for inputting samples - Two barcoded samples (01, 02) corresponding to those in the sample sheets. These data are from the example test dataset for `wf-teloseq`, trimmed down to 1000 reads.

## Test cases
- **defaults**:
    > Test the default command, which is `--fastq HG002_small_test.fastq.gz --reference data/HG002qpMP_reference.fasta.gz`. Then compares the output summary csv to the expected output.
- **bam_in**:
    > Runs the workflow, ingesting the sample bams in test_data/samples.
- **skip_mapping**:
    > Runs the single fastq file with the `--skip_mapping` flag applied.
- **non_telomeric_reads**:
    > Runs the workflow, inputting reads that are non-telomeric. This doesn't fail, it skips every step after `check_reads`. 
- **missing_reference**:
    > Runs the workflow with a missing, reference, expected to fail. Doesn't check for this, so more explodes than exits gracefully.
- **mapq61**:
    > Sets the --mapq flag to 61, which correctly causes the pipeline to fail.
- **sample_sheet_same**
    > Runs with a `--sample-sheet` set, where the reference in the sample sheet for each sample is the same reference.
- **sample_sheet_different**
    > Runs with a `--sample-sheet` set, where the reference in the sample sheet for each sample is different.