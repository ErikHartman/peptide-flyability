import subprocess
import pandas as pd
from Bio import SeqIO

SLICES = 5
records = []
for record in SeqIO.parse(
    "/Users/erikhartman/dev/peptide-flyability/raw_data/UP000005640_9606.fasta",
    format="fasta",
):
    records.append(record)

N = len(records) // SLICES
for i in range(SLICES):
    with open(f"/Users/erikhartman/dev/peptide-flyability/raw_data/UP000005640_9606_slice_{i}.fasta", "w") as f:
        for record in records[i * N : (i + 1) * N]:
            SeqIO.write(record, f, format="fasta")
    cmd = [
        "sage",
        "search.json",
        "-o",
        f"semi_{i}",
        "-f",
        f"/Users/erikhartman/dev/peptide-flyability/raw_data/UP000005640_9606_slice_{i}.fasta",
        "--write-pin",
        "HeLa_chytry_HCD_1.mzML.gz",
    ]
    subprocess.run(cmd)


dfs = []
for i in range(SLICES):
    dfs.append(pd.read_csv(f"semi_{i}/results.sage.pin", sep="\t"))

pd.concat(dfs).sort_values(by="ln(hyperscore)", ascending=False).drop_duplicates(
    subset=["FileName", "ScanNr"], keep="first"
).to_csv("sliced.pin", sep="\t")
