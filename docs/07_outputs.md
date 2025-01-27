Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | wf-teloseq-report.html | Report for all samples | aggregated |
| Tool versions | versions.txt | A CSV with per row tool and version. | aggregated |
| Parameters from workflow | params.json | A json of all parameters selected in workflow. | aggregated |
| Chromosome arm statistics using strict filtering | {{alias}}/results/high_filtered_chr_arm_coverage.csv | Chromosome arm statistics using strict filtering. | per-sample |
| Chromosome arm statistics using lenient filtering | {{alias}}/results/low_filtered_chr_arm_coverage.csv | Chromosome arm statistics using lenient filtering. | per-sample |
| Per telomere read statistics using strict filtering | {{alias}}/results/high_filtered_per_read_telomere_length.csv | Per telomere read statistics using strict filtering. | per-sample |
| Per telomere read statistics using lenient filtering | {{alias}}/results/low_filtered_per_read_telomere_length.csv | Per telomere read statistics using lenient filtering. | per-sample |
| Per telomere read statistics using no filtering | {{alias}}/results/no_filtered_per_read_telomere_length.csv | Per telomere read statistics using no filtering. | per-sample |
| Per telomere read statistics of bulk | {{alias}}/results/sample_raw_per_read_telomere_length.csv | Per telomere read statistics of bulk. | per-sample |
| Trimmed telomere reads | {{alias}}/reads/reads_trimmed.fastq | Adapter trimmed identified telomere reads. | per-sample |
| Read statistics | {{alias}}/results/output.csv | Read stats for each filter. | per-sample |
| Aligned reads strict | {{alias}}/alignments/{{alias}}.tagged.bam | Aligned reads after strict filtering. | per-sample |
| De novo contig naming | {{alias}}/reference_naming/naming_summary.csv | Summary of blast results used for contig naming. | per-sample |
| Reference used for alignment | {{alias}}/alignments/filtered_reference.fasta | Reference used for alignment after extraction. | per-sample |
