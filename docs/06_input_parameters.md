### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| doublestranded | boolean | double stranded protocol | If selected then the telomere reads are identified on both strands, then the second strand is reverse complemented so that all reads are orientated telomere first and pipeline continues as normal but with the other strand reads used. | False |
| reference | string | Reference genome of the sequenced sample. | Reference genome of the sequenced material in fasta format. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### TeloSeq Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| skipmapping | boolean | Skip mapping step for just sample only telomere length | If selected then the workflow will not run the mapping step but measure telomere length just on the unmapped telomere identified reads. | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| adapters_to_trim | string | Comma-delimited list of sequencing adapters to be trimmed |  | CACCCTAACCCTAACCCTAACC,CAACCCTAACCCTAACCCTAAC,CCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCT,CTAACCCTAACCCTAACCCTAA,CCCCTAACCCTAACCCTAACCC |
| filter_error_motifs | string | Comma-delimited list of telomere error motifs for filtering | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in a telomere read, it is dropped from the analysis. | GTATAG,CGCGCGCG,CCACCG,AGCGACAG,ATAAGT,CCTCGTCC,TATAGT,AGTACT,GAGTCC,TATAGT,TATACA,TGGTCC,CTCTCCTCT |
| filter_error_motifs_max_count | integer | Number of erroneous telomere k-mers to remove read. | If more than this number of telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in the telomeric region of a read, it is dropped from the analysis. | 5 |
| filter_error_motifs_window_size | integer | Size of window for filtering based on telomere error motifs. | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of this length in the telomeric region of a read, it is dropped from the analysis. | 500 |
| mapq | integer | mapping quality filter parameter | Mapping quality used to filter the bam file | 4 |
| read_quality | integer | read quality filter parameter | Read quality used to filter raw reads | 9 |
| restriction_site | string | enzyme cut site | Restriction enzyme cut site used for filtering reads that are not close to this site for the strict setting | GATATC |
| cov_4cluster | integer | minimum read number in clustering to produce contig, recommend 30 for 220x and 20k telomere reads and 8 for 3k telomere reads. | The minimum read number to be used for filtering after using the clustering algorithm | 7 |
| mincoverage | integer | minimum read number for coverage of reference contig, default is 20% of telomere read average for 92 chr arms | The minimum telomere coverage of chromosome arms to be taken to the final results and plots is calculated as 20% of the average chr arm coverage but can be overridden by giving a value here | -1 |


