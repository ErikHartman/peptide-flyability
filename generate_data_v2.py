import pandas as pd
import numpy as np
import re
from itertools import combinations
from dpks import QuantMatrix

class FeatureGenerator:
    """
    Positional encoding
    
    """

    def __init__(self, remove_unimod : bool = False) -> None:
        self.remove_unimod = remove_unimod
        self.datasets = []

    def add_dataset(self, data: pd.DataFrame, design: pd.DataFrame):
        data, design = self.preprocess(data, design, self.remove_unimod)
        self.datasets.append((data, design))

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

                        a_feat = create_positional_matrix(a)
                        b_feat = create_positional_matrix(b)

                        a_int = protein_data.set_index("PeptideSequence").loc[a, sample]
                        b_int = protein_data.set_index("PeptideSequence").loc[b, sample]

                        feature_vector = flatten_and_concatenate([a_feat, b_feat])
                        int_diff = np.sum(np.log2(a_int)) / np.sum(np.log2(b_int))

                        feat_id = f"{sample}_{protein}_{a}_{b}"

                        X[feat_id] = feature_vector
                        y[feat_id] = int_diff

                break # Here for speed. Only 1 sample for now.

            Xs.append(pd.DataFrame(X, index=get_index()).T)
            ys.append(pd.Series(y))
            dataset_index += 1

        print("Concatenating")
        X = pd.concat(Xs)
        y = pd.concat(ys)
        return X, y

    

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


def normalize_sequence(sequence):
    max_length = 50
    return sequence.ljust(max_length, 'X')[:max_length]


def create_positional_matrix(sequence):
    sequence = remove_unimod(sequence)
    normalized_seq = normalize_sequence(sequence)
    matrix = np.zeros((len(normalized_seq), len(amino_acids)))
    for pos, aa in enumerate(normalized_seq):
        matrix[pos, amino_acid_index[aa]] += 1
    return matrix


def flatten_and_concatenate(matrices):
    vectors = [matrix.flatten() for matrix in matrices]
    return np.hstack(vectors)

def get_index():
    index = []
    for pos in range(50):
        for aa in amino_acids:
            index.append(f"{pos}_{aa}")
    return [i + "_1" for i in index] + [i + "_2" for i in index] 

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
    "X"
]
amino_acid_index = {aa: idx for idx, aa in enumerate(amino_acids)}
