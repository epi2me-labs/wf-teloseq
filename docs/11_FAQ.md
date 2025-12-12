<!---Frequently asked questions, pose any known limitations as FAQ's.--->

- _Can I use TeloSeq and `wf-teloseq` on non human samples?_ Yes, although Telo-Seq has been tested on human data, it should work well on species with similar telomere repeat sequences.

- _How can I prepare my own reference for alignment?_ The workflow makes several assumptions regarding the sequences present in the provided reference database (FASTA file). For preparing a reference for use with `wf-teloseq`, the workflow requires:
    - Telomeres on opposing arms of a chromosome MUST be provided as *separate* reference sequences (named as below). It is not permissible to provide a reference sequence of a complete, end-to-end, chromosome.
    - Reference contigs should be trimmed to the first ecoRV cutsite after the end of the telomeric repeat region. The cutsite sequence is `GATATC`. 
    - Reference contig sequences must be orientated 3&#8242; -> 5&#8242;.  
      This will place the telomere on the left of the contig, followed by subtelomere.
    - Reference sequences must be named in the form `{chromosome}_{mat/pat}_{p/q}` (or `{chromosome}_{p/q}`, see below).
      For example, `chr1_mat_q`, meaning the telomere of the q-arm of chromosome 1.
      Note: The workflow explicitly checks the presence of `mat` or `pat` in the contig name to assign haplotype.
    - If a diploid reference is not available, or results at this level of detail not required, the haplotype identifiers can be omitted. For example, `chr1_q`.

- _What happened to `Sample_raw_Per_Read_telomere_length.csv` from the beta version of the workflow?_ This file was removed, as all per-read information is contained in the output BAM file. Users are expected to use tools such as [pysam](https://pysam.readthedocs.io/en/stable/) if they wish to perform their own per-read analyses. Alternatively raise a feature request if there's a feature you'd like to see added to the workflow reports.

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-teloseq/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 