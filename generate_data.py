import pandas as pd
import numpy as np
import re
from itertools import combinations
from dpks import QuantMatrix


from Bio.SeqUtils.ProtParam import ProteinAnalysis


class FeatureGenerator:

    def __init__(self, remove_unimod : bool = False) -> None:
        self.index_no_unimod = ["gravy", "mw", "arom", "ip", "length"] + list(
            amino_acids
        )
        self.remove_unimod = remove_unimod
        self.datasets = []

    def add_dataset(self, data: pd.DataFrame, design: pd.DataFrame):
        data, design = self.preprocess(data, design, self.remove_unimod)
        self.datasets.append((data, design))
        self.unique_unimods = get_unique_unimods(data, "PeptideSequence")
        index = self.index_no_unimod + [unimod for unimod in self.unique_unimods]
        self.total_index = [f"{i}_1" for i in index] + [f"{i}_2" for i in index]
        print(self.unique_unimods)

    def generate_pair_matrix(self, topn: int = None):
        Xs = []
        ys = []
        dataset_index = 0
        for data, design in self.datasets:

            data = data.copy()
            design = design.copy()

            print(f"Dataset {dataset_index}")

            data: pd.DataFrame
            design: pd.DataFrame
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
                        protein_data = protein_data.sort_values(
                            sample, ascending=False
                        )[0:topn]

                    peptides = protein_data["PeptideSequence"].tolist()

                    peptide_pairs = get_pairs(peptides)

                    for ab in peptide_pairs:

                        a, b = ab

                        a_feat = self.generate_peptide_feature_vector(a)
                        b_feat = self.generate_peptide_feature_vector(b)

                        a_int = protein_data.set_index("PeptideSequence").loc[a, sample]
                        b_int = protein_data.set_index("PeptideSequence").loc[b, sample]

                        feat_diff = np.array(a_feat + b_feat)
                        int_diff = np.log2(a_int) / np.log2(b_int)

                        feat_id = f"{sample}_{protein}_{a}_{b}"

                        X[feat_id] = feat_diff
                        y[feat_id] = int_diff
                break # Here for speed. Only 1 sample for now.
            Xs.append(pd.DataFrame(X, index=self.total_index).T)
            ys.append(pd.Series(y))
            dataset_index += 1
        print("Concatenating")
        X = pd.concat(Xs)
        y = pd.concat(ys)
        return X, y

    def generate_peptide_feature_vector(self, peptide: str):
        peptide_no_unimod = remove_unimod(peptide)
        analyzed_seq = ProteinAnalysis(peptide_no_unimod)
        gravy = analyzed_seq.gravy()
        molecular_weight = analyzed_seq.molecular_weight()
        amino_acids_percent = analyzed_seq.get_amino_acids_percent()
        aromaticity = analyzed_seq.aromaticity()
        isoelectric_point = analyzed_seq.isoelectric_point()
        length = len(peptide_no_unimod)
        feature_vector = [
            gravy,
            molecular_weight,
            aromaticity,
            isoelectric_point,
            length,
        ]
        feature_vector += [amino_acids_percent[aa] for aa in amino_acids]
        feature_vector += self.count_unimods(peptide)
        return feature_vector

    def preprocess(self, data, design, remove_mod=True):
        data["PeptideSequence"] = data["PeptideSequence"] + ";" + data["Protein"]
        data = (
            QuantMatrix(data, design)
            .quantify(
                method="top_n",
                top_n=6,
                summarization_method="sum",
                level="PeptideSequence",
            )
            .to_df()
            .replace(0, np.nan)
        )
        data[["PeptideSequence", "Protein"]] = data["PeptideSequence"].str.split(
            ";", expand=True
        )
        if remove_mod:
            data["PeptideSequence"] = data["PeptideSequence"].str.replace(")", ">>")
            data = data[~data["PeptideSequence"].str.contains(">>")]
        return data, design

    def count_unimods(self, peptide: str):
        unimods = re.findall(r"\(UniMod:\d+\)", peptide)
        unimod_counts = {u:0 for u in self.unique_unimods}
        for unimod in unimods:
            unimod_counts[unimod] += 1

        return list(unimod_counts.values())

def get_unique_unimods(df: pd.DataFrame, column_name: str):
    unique_unimods = set()
    unimod_pattern = re.compile(r"\(UniMod:\d+\)")

    for peptide in df[column_name]:
        unimods = unimod_pattern.findall(peptide)
        unique_unimods.update(unimods)

    return list(unique_unimods)





def remove_unimod(peptide_sequence: str):
    peptide_sequence = re.sub(r"\(UniMod:\d+\)", "", peptide_sequence)
    peptide_sequence = re.sub(r"\(\d+\)", "", peptide_sequence)
    return peptide_sequence


def get_pairs(objects: list):
    return list(combinations(objects, 2))


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
