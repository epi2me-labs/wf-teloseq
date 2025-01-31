# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]
### Added
- Added coefficient of variance into report.
- Added naming of contigs using EcoRV pangenome for Human but user can supply their own pangenome. Limited to those can name with certainty, which is approximately halve of paternal/maternal arms due to similarity that also makes mapping to a pangenome problematic.
- Addition of reference to meta data allows for use of different references per sample via --sample_sheet option and use of the same keys throughout pipeline.

### Changed
- Re-written the de novo guided route to be fully de novo.
- Lenient filter extended to 2000 bp and if cut site before this then switch to cut site. This reduces some very low level mismapping tolerance. 
- Lenient and strict filter changed to high and low stringency renaming.

### Removed
- Removed the "de novo" guided route. This will be added back in v0.1. 


## [v0.0.4]
### Changed
- Bug fix for multi-sample folder input where samples results where getting mixed.


## [v0.0.3]
### Added
- Paramater for min coverage if don't want 20% of average telomere read number over 92 chromosomes (default)
- Paramater if you have protocol to capture both strands
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
