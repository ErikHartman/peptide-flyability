import pandas as pd
import numpy as np
import re
from itertools import combinations



from Bio.SeqUtils.ProtParam import ProteinAnalysis


amino_acids = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]

index = ["rt", "gravy", "mw", "arom", "ip", "charge", "length"] + list(amino_acids)


def remove_unimod(peptide_sequence: str):
    peptide_sequence = re.sub(r"\(UniMod:\d+\)", "", peptide_sequence)
    peptide_sequence = re.sub(r"\(\d+\)", "", peptide_sequence)
    return peptide_sequence


def generate_pair_matrix(data: pd.DataFrame, design: pd.DataFrame, topn : int = None):

    samples = design["sample"]
    proteins = data["Protein"].unique()

    X = {}
    y = {}

    for sample in samples:
        for protein in proteins:
            protein_data = data.loc[data["Protein"] == protein].copy()

            protein_data: pd.DataFrame
            protein_data = protein_data.dropna(subset=sample)
            if len(protein_data) == 0:
                continue

            if topn:
                protein_data = protein_data.sort_values(sample, ascending=False)

            peptides = protein_data["Precursor"].tolist()

            rts = protein_data["RetentionTime"].tolist()
            rts = {peptide: rt for peptide, rt in zip(peptides, rts)}

            charges = protein_data["Charge"].tolist()
            charges = {peptide: charge for peptide, charge in zip(peptides, charges)}

            peptide_pairs = get_pairs(peptides)

            for ab in peptide_pairs:

                a, b = ab

                a_feat = generate_peptide_feature_vector(a, rts[a], charges[a])
                b_feat = generate_peptide_feature_vector(b, rts[b], charges[b])

                a_int = protein_data.set_index("Precursor").loc[a, sample]
                b_int = protein_data.set_index("Precursor").loc[b, sample]

                feat_diff = np.array(a_feat) - np.array(b_feat)
                int_diff = np.log2(a_int) / np.log2(b_int)

                feat_id = f"{sample}_{protein}_{a}_{b}"

                X[feat_id] = feat_diff
                y[feat_id] = int_diff

        return pd.DataFrame(X, index=index).T, pd.Series(y)


def generate_peptide_feature_vector(peptide: str, rt: float, charge: float):
    peptide = remove_unimod(peptide)
    analyzed_seq = ProteinAnalysis(peptide)
    gravy = analyzed_seq.gravy()
    molecular_weight = analyzed_seq.molecular_weight()
    amino_acids_percent = analyzed_seq.get_amino_acids_percent()
    aromaticity = analyzed_seq.aromaticity()
    isoelectric_point = analyzed_seq.isoelectric_point()
    length = len(peptide)
    feature_vector = [
        rt,
        gravy,
        molecular_weight,
        aromaticity,
        isoelectric_point,
        charge,
        length,
    ]
    feature_vector += [amino_acids_percent[aa] for aa in amino_acids]
    return feature_vector


def get_pairs(objects : list):
    return list(combinations(objects, 2))