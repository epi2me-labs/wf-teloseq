# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- Addition of reference to meta data allows for use of different references per sample via --sample_sheet option and use of the same keys throughout pipeline.
- New tags for telomere boundary, quality control outcomes and haplotype for each read.
- Quality control outcomes are listed in the workflow documentation. Added at this stage are:
    - Too Short - The read is less than two window widths long
    - Too Few Repeats - Read did not meet minimum number of telomere repeats across whole read 
    - Start Not Repeats - The start of the read is not comprised of a sufficient number of repeats
    - Too Close End - Telomere boundary is half window width from end of read
    - Too Close Start - The boundary is within one window width of the start of the read
    - Low Sub-Telo Qual - The median quality of the basecalled read after the telomere boundary is too low.
    - Too Errorful - Too many known basecalling error motifs are present in the telomere.
    - Bad Alignment - The read has a low gap compressed identity to the reference, or a too low mapping quality.
    - Telomere Only - The composition of the post telomere boundary sequence for non-pverlapping `CCC` is more than 25%.
- `--alignment_threads` parameter to control the number of threads used by minimap2. Defaults to 6. 

### Changed
- The filtering of the input data has been refactored, and boundary detection has been altered. This has removed several processes, collapsing them all down into one.
- Added custom script for boundary detection, which detects the coordinate of large shifts in the density of telomeric repeats.
- The alignment process has had it's maximum memory directive raised to 7Gb from 2Gb. We were seeing some occurrences of the process being killed for exceeding the memory cap on larger datasets when using higher thread counts.
- Documentation has been updated and rewritten.
- Output BAM now contains all input data, even unmapped and those lacking detectable telomeric boundaries, allowing for further analysis by the user.
- Reads which fail QC will still have any detected telomere boundary tagged on the record, in order to enable further downstream investigation.
- Updated to wf-template v5.6.0, changing:
    - Reduce verbosity of debug logging from fastcat which can occasionally occlude errors found in FASTQ files during ingress.
    - Log banner art to say "EPI2ME" instead of "EPI2ME Labs" to match current branding. This has no effect on the workflow outputs.

### Removed
- Removed the "de novo" guided route.
- Removed the Greider telomere boundary detection code.

### Fixed
- Pipeline now correctly handles empty input files in a multi sample run. Samples with no input will be filtered out, and a warning will be logging per sample removed.
- Samples with no good telomere data will be carried through the workflow, and correctly tagged. The report will contain basic statistics about the sample, but all telomere related stats and plots will contain explainers that there was no data.
- Minimum and Maximum values displayed in table are now the same as those used in the Boxplots (within 1.5 x interquartile range of the median).
- Updated to wf-template v5.6.0, fixing:
    - dacite.exceptions.WrongTypeError during report generation when barcode is null.
    - Sequence summary read length N50 incorrectly displayed minimum read length, it now correctly shows the N50.
    - Sequence summary component alignment and coverage plots failed to plot under some conditions.

## [v0.0.4]
### Changed
- Bug fix for multi-sample folder input where samples results where getting mixed.


## [v0.0.3]
### Added
- Parameter for min coverage if don't want 20% of average telomere read number over 92 chromosomes (default)
- Parameter if you have protocol to capture both strands
- Adapter trimming (note not precise in terms of sequence removes as may take 1-2 telomere bases in removal but removes all adapter presence up to telomere)
- Parameter for enzyme cut site sequence

### Changed
- Adopted telomere boundary code adapted from Greider lab published code (MIT licence). Both telomere boundary sites are reported including prior method in output but Greider taken for plots.

## [v0.0.2]
### Added
- Functionality for mapping and de novo guided reference.
- Extra pathway for when performed denovo reference-guided route and then want to add manual curated contigs before mapping.
- Added read quality filter set default to 0 but can be set to recommended 9.
- Put in minimum coverage threshold for reporting chr arm coverage and reference contig creation using 25% of telomere reads identified over 92 chr.

### Changed
- Refactored to Epi2me standards.
- Report rejigged to handle multiple samples better and also output general read stats.
- Output structure changed to be within the output folder structure required by Epi2me.
- The documentation has been updated.
- Changed lenient filter so if subtelomere is very short that don't just remove all those reads, added catch if last telomere plus 80 is beyond contig length.
- Reduced mapq from 10 to 4, 10 was too restrictive.
- Changed strict filter to remove length threshold and just add two references by name that have high repetitive content.

## [v0.0.1]
- Initial commit
