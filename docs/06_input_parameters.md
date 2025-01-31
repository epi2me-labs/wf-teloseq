### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| reference | string | Reference genome of the sequenced sample. |  |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .FASTQ(.gz) file or directory of .FASTQ(.gz) files. |  |  |


### TeloSeq Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| skip_mapping | boolean | Skip mapping step for just sample only telomere length | If selected then the workflow will not run the mapping step but measure telomere length just on the unmapped telomere-subtelomere identified reads. | False |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| adapters_to_trim | string | Comma-delimited list of sequencing adapters to be trimmed up to but not including. | Uses telomere so all adapter is trimmed off the ends of the telomere containing end of reads. | CACCCTAACCCTAACCCTAACC,CAACCCTAACCCTAACCCTAAC,CCTAACCCTAACCCTAACCCTA,CCCTAACCCTAACCCTAACCCT,CTAACCCTAACCCTAACCCTAA,CCCCTAACCCTAACCCTAACCC |
| filter_error_motifs | string | Comma-delimited list of telomere error motifs for filtering | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in a telomere read, it is dropped from the analysis. | GTATAG,CGCGCGCG,CCACCG,AGCGACAG,ATAAGT,CCTCGTCC,TATAGT,AGTACT,GAGTCC,TATAGT,TATACA,TGGTCC,CTCTCCTCT |
| filter_error_motifs_max_count | integer | Number of erroneous telomere k-mers to remove read. | If more than this number of telomere error motifs are found within a window of size `--filter_error_motifs_window_size` in the telomeric region of a read, it is dropped from the analysis. | 5 |
| filter_error_motifs_window_size | integer | Size of window for filtering based on telomere error motifs. | If more than `--filter_error_motifs_max_count` telomere error motifs are found within a window of this length in the telomeric region of a read, it is dropped from the analysis. | 500 |
| mapq | integer | Mapping quality filter parameter | Mapping quality used to filter the BAM file | 4 |
| min_coverage_percent | integer | Minimum percentage coverage of total telomere reads per chr arm for filtering | Used in minimum coverage calculation for telomere reads using number of chromosome arms | 15 |
| min_length | integer | Minimum read length for filtering | Used in initial filtering of reads. | 100 |
| exclude_chr_from_naming | string | Exclusion naming list | Used in chromosome naming to exclude from naming as too similar to others | chr9q chr11p chr13p chr15p chr19p chr6p chr20q chr3q chr6q chr7p chr14p |
| motif_threshold | integer | Number of motif of repeat to identify chr21p | Chr21p can be identified by the repeat composition rather than blastn | 40 |
| max_sstart | integer | Max starting location of blastn hit | This ensures blasthits are not part of the reference but spans the beginning | 20 |
| min_pident | number | Minimum percent identity from the blastn | Used in the de novo naming to filter poor hits | 95 |
| exp_chr_num | integer | Expected number of chr | Used in minimum coverage calculation for telomere reads with percentage minimum covereage number | 92 |
| naming_file | string | pangenome used to name de novo contigs |  |  |
| read_quality | integer | Read quality filter parameter |  | 9 |
| restriction_site | string | Enzyme cut site | Restriction enzyme cut site used for filtering reads that are not close to this site for the strict setting | GATATC |
| beyond_cut | integer | Amount of reference to include beyond cut site for each contig |  | 300 |
| telomere_extension | integer | Addition of telomere to reference avoid mismapping | Reads primary and secondary alignments can be incorrect if one similar reference contig has longer telomere that allows read to match to it better than the other contig | 4000 |
| telomere_margin | integer | Distance from telomere boundary to use to remove reads whose 3' end do not map up to this chromosome position. |  | 2000 |
| min_coverage | integer | minimum read number for coverage of reference contig, default is 20% of telomere read average for 92 chr arms | The minimum telomere coverage of chromosome arms to be taken to the final results and plots is calculated as 20% of the average chr arm coverage but can be overridden by giving a value here | -1 |


