### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| reference | string | Reference genome of the sequenced sample, if not specified a human derived telomere reference set will be used. |  |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .FASTQ(.gz) file or directory of .FASTQ(.gz) files. |  |  |


### TeloSeq Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| skip_mapping | boolean | Perform alignment to assign haplotypes to telomeric reads. | Use `--skip_mapping` if there is no suitable reference available. Only a bulk estimate of telomere length will be calculated. | False |
| alignment_threads | integer | Set max number of threads to use for alignment. |  | 8 |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| mapq | integer | Mapping quality filter parameter | Alignments to the reference with a mapping quality (MAPQ) score lower than this value will be discarded from further processing. | 4 |
| min_length | integer | Minimum read length for filtering | Used in initial filtering of reads. This removes all reads which are either noise, artifacts or partial, and which cannot be full telomeric sequences. Setting this number too high runs the risk of artificially inflating calculated telomere lengths. | 100 |
| read_quality | integer | Reads with a mean quality score lower than this value will be filtered out prior to analysis. |  | 9 |
| restriction_site | string | Enzyme cut site | Restriction enzyme cut site used for filtering reads that are not close to this site for the strict setting | GATATC |
| beyond_cut | integer | Amount of reference to include beyond cut site for each contig |  | 300 |
| telomere_extension | integer | Addition of telomere to reference avoid mismapping | Reads primary and secondary alignments can be incorrect if one similar reference contig has longer telomere that allows read to match to it better than the other contig | 4000 |
| telomere_margin | integer | Distance from telomere boundary to use to remove reads whose 3' end do not map up to this chromosome position. |  | 2000 |


