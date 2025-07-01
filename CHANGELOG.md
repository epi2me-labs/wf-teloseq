# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased
### Changed
- Updated to wf-template v5.6.2 to maintain compliance with our latest wf-template standard: this does not impact the workflow.

### Fixed
- Sample plots were previously linked to incorrect barcodes in the 'Read summary plots' dropdown selector of the report. This issue has been fixed, and the plots are now correctly associated with their respective sample aliases. This fix does not alter report plotting, only the sample dropdown ordering.
- KDE calculations previously crashed the workflow when only a single read passed filtering. This issue has now been handled, and KDE plots will only display when two or more reads pass filtering and have different lengths.

## [v1.0.0]

wf-teloseq is a new workflow designed to analyse the output of the Telo-Seq protocol.
Our v1.0.0 release supersedes the experimental beta release and is focused on a speedy analysis of telomeric lengths, either using a sample matched reference to measure telomere lengths per chromosome arm, or, if such a reference is not available, providing a per sample length profile for all input sequences which were identified as telomeric.
We look forward to seeing how this supports your research and discovery using the Telo-Seq protocol.
