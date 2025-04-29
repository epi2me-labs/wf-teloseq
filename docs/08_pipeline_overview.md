The workflow is composed of two steps. 
Firstly, reads undergo initial basic length and quality control filtering, followed by analysis to determine the telomeric boundary, and filtering to remove reads where the boundary is likely incorrect.
If alignment is enabled (default), all reads (irrespective of boundary determination) are then aligned to a reference. Reads with primary alignments with good telomere boundaries are then aggregated, before statistics about the estimated lengths of each contigs telomeres are generated. 

## 1. Input and Sequence preparation.

`wf-teloseq` requires basecalled nanopore reads in FASTQ or BAM format as input.
We recommend basecalling with super accuracy (SUP) models.
High accuracy (HAC) basecalling will also work but will result in a slightly reduced number of telomeric reads of sufficient quality for analysis.
For best results, we recommend a minimum of 1000 telomere reads per sample for alignment based analyses.
If a reference is not provided via the `--reference` or `--sample_sheet` parameter, then a telomere reference constructed from HG002 (`data/HG002qpMP_reference.fasta.gz`) is used by default for alignment.

See our [Know-How document](https://community.nanoporetech.com/knowledge/know-how/TELO-seq) for further information.


## 2. Read filtering and tagging
Filters are listed in the order that they are applied.
The first two filters are tunable with the params `--min_length` and `--read_quality`. 
Other filters are currently not tunable via parameters.

Reads which fail a check are tagged with a `qc:Z:<Tag>` tag, in either the output FASTQ if not performing alignment or the output BAM.
Tags for corresponding failure modes are listed below.
Reads which pass all filtering checks are tagged with `qc:Z:Good`.

| Filter               | Description                                                                 | Tag               |
|----------------------|---------------------------------------------------------------------------|------------------|
| **Minimum read length** | Reads under 100 bases long are discarded. Removed entirely.             | N/A              |
| **Minimum Q score**  | Reads with a Mean Q score <9 are discarded. Removed entirely.            | N/A              |
| **Too Short**        | The read was too short (less than 160 bases). Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooShort        |
| **Too Few Repeats**  | There were fewer than 20 telomeric repeat motifs across the entire read. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooFewRepeats   |
| **Too Close Start**    | The telomeric boundary is too close (within 60 bases) of the Start of the read. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooCloseStart     |
| **Too Close End**    | The telomeric boundary is too close (within 60 bases) of the end of the read. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooCloseEnd     |
| **Start Not Repeats** | The first 30% of the read is not 80% repeats. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | StartNotRepeats |
| **Low Sub-Telo Qual** | The mean basecall Q score of the region after the boundary is below a default value of 9. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | LowSubTeloQual  |
| **Too Errorful**     | At least 5 known basecall error motifs have been observed within a 500 base pair frame in the subtelomere. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | TooErrorful |
| **Telomere only**     | Sequence after telomere boundary is CCC rich, and so is most likely still actually telomeric repeat sequence, meaning the detected boundary is unreliable. | TelomereOnly |

The following filter is only applied if alignment is performed:
| Filter               | Description                                                                 | Tag               |
|----------------------|---------------------------------------------------------------------------|------------------|
| **Bad Alignment**    | Query read has a low Gap compressed identity to the reference (\< 0.8), or a low (\< 20) mapping quality. Read is excluded from further analysis, is tagged as failing QC, remaining in output BAM. | BadAlign        |

All reads have a `tl:i` tag set, which represents the estimated telomere length. 
For reads which a boundary CANNOT be determined, the value of this tag is set to -1.

It is after this stage that a tagged FASTQ file is available per sample, with the tags present in the header.

## 3. Alignment.
By default, ALL reads are aligned by minimap2 (excepting reads which fail the initial minimum read length and Q score filters) after initial processing. 
Alignment can be skipped by setting `--mapping false` when running the pipeline.

### Default reference genome

The telomere reference provided with the pipeline is based on the [HG002 telomere-to-telomere reference genome](https://github.com/marbl/hg002).
Note that in this reference, the sub-telomeric sequence up to the EcoRV cut site in chromosome 13 (paternal) P arm is identical to chromosome 22 (paternal) P arm, so there is only one representative sequence present in the file, with a sequence name of `chr13_22_PATERNAL_P`.
The pipeline does not separate out these two arms based upon telomere length from the single contig provided, which means the estimated telomere length assigned to this contig is likely inaccurate.

### Processing of alignments
After alignment, the alignments are used to aggregate stats based on the reference contigs.
One final filtering step is first performed here, which is based on the Gap compressed identity between the query sequence and the target reference, as calculated by minimap2.
Any read with an identity of less than 0.8 and a mapping quality of less than 20 is excluded from the estimated telomere length statistics.
Only primary alignments are used in this aggregation, and only those from reads which have the `qc:Z` tag set with a value of `Good`, indicating they passed all initial filtering steps.

After this stage a tagged BAM file containing both mapped and unmapped reads is produced and output per sample.
