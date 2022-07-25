# This script was used to combine the exon junctions with the half-life datasets

import pandas as pd


ej = pd.read_csv("D:/Documents/Bioinformatik/Master/ML/project_1/data/GTEx_exon_junctions.txt", delimiter="\t")
ej.drop(["UTR5_len"], axis=1, inplace=True)
ej["Exon_Junctions_In_Full_Sequence"] = ej["Exon_Junctions_In_Full_Sequence"].str.replace(',', ';')
ej.fillna('', inplace=True)

# hl = pd.read_csv("D:/Documents/Bioinformatik/Master/ML/project_1/data/tissue_hl.csv", delimiter=",")
hl = pd.read_csv("D:/Documents/Bioinformatik/Master/ML/project_1/data/genomic_sequence_plus_features_hl_all_tissues.csv", delimiter=",")
hl.columns.values[0] = "GeneID"
hl["GeneID"] = hl["GeneID"].str[:15]

combined = pd.merge(ej, hl, on="GeneID")
combined.rename(columns={'Exon_Junctions_In_Full_Sequence': 'Exon_Junctions'}, inplace=True)

combined.to_csv("D:/Documents/Bioinformatik/Master/ML/project_1/data/genomic_sequence_plus_features_hl_all_tissues_with_ss.csv", index=False)
