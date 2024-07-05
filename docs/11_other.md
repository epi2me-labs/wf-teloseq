## Additional

Telo-Seq has only been tested on human data, but we expect it to work well on species with similar telomere repeat sequences.

Telomere sequences are highly repetitive. Sometimes basecalling errors cause dense mis-basecalled repeats to occur in Telo-seq reads. These get soft-clipped off by the aligner in downstream analysis, which could result in underestimation of telomere lengths. `wf-teloseq` identifies such reads and excludes them from telomere length estimation.


# Acknowledgements

This project uses code in the "telomerewindowV1.py" script from the following source:

- **Original Author:** Ramin Kahidi
- **Original Repository:** https://github.com/GreiderLab/TeloBP

The code used is licensed under the MIT License, which can be found in the original repository and the header of the script.

Publication: Karimian K, Groot A, Huso V, Kahidi R, Tan KT, Sholes S, Keener R, McDyer JF, Alder JK, Li H, Rechtsteiner A, Greider CW. Human telomere length is chromosome specific and conserved across individuals. 2024 Jan 13:2023.12.21.572870. doi: 10.1101/2023.12.21.572870. Update in: Science. 2024 May 3;384(6695):533-539. doi: 10.1126/science.ado0431. PMID: 38187739; PMCID: PMC10769321.


# Epi2Me

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.