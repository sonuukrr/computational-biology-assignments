from pathlib import Path

import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
ANNOTATION_FILE = BASE_DIR / "human_gene_annotation.tsv.gz"
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)


def normalize_chromosome(chromosome):
    chromosome = str(chromosome)
    return chromosome if chromosome.startswith("chr") else f"chr{chromosome}"


def main():
    annotation = pd.read_csv(
        ANNOTATION_FILE,
        sep="\t",
        compression="gzip",
        dtype={"chromosome_name": str},
        low_memory=False,
    )
    annotation = annotation.dropna(subset=["chromosome_name", "external_gene_name", "transcription_start_site"])
    annotation["chrom"] = annotation["chromosome_name"].map(normalize_chromosome)
    annotation["tss"] = annotation["transcription_start_site"].astype(int)
    annotation["strand_symbol"] = annotation["strand"].map({1: "+", -1: "-"})
    annotation = annotation.dropna(subset=["strand_symbol"])

    bed = pd.DataFrame(
        {
            "chrom": annotation["chrom"],
            "start": annotation["tss"].clip(lower=0),
            "end": annotation["tss"].clip(lower=0) + 1,
            "name": annotation.apply(
                lambda row: f"{row['chrom']}@{row['tss']}-{row['tss'] + 1}|{row['external_gene_name']}",
                axis=1,
            ),
            "score": ".",
            "strand": annotation["strand_symbol"],
        }
    ).drop_duplicates()

    promoter_start = []
    promoter_end = []
    for row in bed.itertuples(index=False):
        if row.strand == "+":
            promoter_start.append(max(0, row.start - 500))
            promoter_end.append(row.end)
        else:
            promoter_start.append(row.start)
            promoter_end.append(row.end + 500)

    promoters = bed.copy()
    promoters["start"] = promoter_start
    promoters["end"] = promoter_end

    bed_path = RESULTS_DIR / "human_tss.bed"
    promoter_path = RESULTS_DIR / "human_tss_500bp_promoters.bed"
    bed.to_csv(bed_path, sep="\t", index=False, header=False)
    promoters.to_csv(promoter_path, sep="\t", index=False, header=False)

    print(f"Wrote {len(bed)} TSS entries to {bed_path}")
    print(f"Wrote strand-aware 500 bp promoter intervals to {promoter_path}")


if __name__ == "__main__":
    main()
