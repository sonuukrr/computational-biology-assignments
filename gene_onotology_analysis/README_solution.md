# Gene Ontology / Promoter Preparation Solution

This folder contains `prepare_promoters.py`, which completes the part that is possible from the provided file `human_gene_annotation.tsv.gz`.

It creates:

- `results/human_tss.bed`: one-base TSS BED records.
- `results/human_tss_500bp_promoters.bed`: strand-aware 500 bp promoter intervals.

Run:

```bash
python prepare_promoters.py
```

The remaining motif search and GO enrichment steps need external files/tools that are not included in this repository:

1. Download `hg38.fa.gz` from UCSC.
2. Install `bedtools`, `emboss`, and Bioconductor `clusterProfiler`.
3. Extract promoter FASTA sequences from `results/human_tss_500bp_promoters.bed`.
4. Search the motif `GCGC..GCGC` using EMBOSS `dreg`.
5. Convert matched promoter names into a gene list.
6. Use `clusterProfiler` to perform GO enrichment on that gene list.
