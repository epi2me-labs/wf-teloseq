## Additional

Telo-Seq has only been tested on human data, but we expect it to work well on species with similar telomere repeat sequences.

Telomere sequences are highly repetitive. Sometimes basecalling errors cause dense mis-basecalled repeats to occur in Telo-seq reads. These get soft-clipped off by the aligner in downstream analysis, which could result in underestimation of telomere lengths. `wf-teloseq` identifies such reads and excludes them from telomere length estimation.

+ [Importing third-party workflows into EPI2ME Labs](https://labs.epi2me.io/nexflow-for-epi2melabs/)

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.