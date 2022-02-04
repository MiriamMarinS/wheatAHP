from Source.utils import comparison, AHP_preparation
import itertools
import ahpy
import pandas as pd
import numpy as np

def criteria1(epitopes_genes_comparisons):
    ahpy_criteria1 = ahpy.Compare(name="Epitopes_genes", comparisons=epitopes_genes_comparisons, precision=3)
    return(ahpy_criteria1)

def criteria2epitope(oligoEpitopes, Epsilon):
    # List of epitopes present in NGS amplicons of alpha- and gamma-gliadins.
    epitopesAlpha = ["p31_L", "DQ2.5_glia_a1a", "DQ2.5_glia_a1b", "DQ2.5_glia_a2", "DQ2.5_glia_a3"]
    epitopesGamma = ["DQ2.5_glia_g2a", "DQ2.5_glia_g4e", "DQ2.5_glia_g4d", "DQ2.5_glia_g4b", "DQ2.5_glia_g1",
                    "DQ2.5_glia_g5", "DQ2.5_glia_g4c", "DQ2.5_glia_g3", "DQ2.5_glia_g4a", "DQ2.5_glia_g2b", "DQ2.5_glia_g4b", "DQ2.5_glia_g4c"]
    epitopesOther = ["DQ2.5_glia_w1", "DQ2.5_hor_2"]
    listEpitopes = epitopesAlpha + epitopesGamma + epitopesOther
    epitopes_with_score = []
    for epitopeID in list(oligoEpitopes.keys()):
        if epitopeID in listEpitopes:
            epitopes_with_score.append(epitopeID)
    comparisons_list_epitopes = list(itertools.combinations(epitopes_with_score, 2)) # All possible combinations of epitopes for pair-wise comparisons.
    dict_criteria2_epitopes = comparison(comparisons_list_epitopes, Epsilon, "epitopes", oligoEpitopes)
    ahpy_criteria2_epitopes = ahpy.Compare(name="Epitopes", comparisons=dict_criteria2_epitopes, precision=3)
    return(ahpy_criteria2_epitopes, epitopes_with_score)

def criteria2gene(rawdata, Epsilon):
    # Independently of the type of amplicon (alpha or gamma), we characterize the amplicons based on the number of epitopes found in them.
    # Amplicon from 0 to 14 epitopes.
    alpha_genes = list(rawdata.iloc[:,np.r_[0,103:110]].columns[1:])
    gamma_genes = list(rawdata.iloc[:,np.r_[0,111:123]].columns[1:])
    list_genes = []
    for number in range(0, 14 + 1):
        alpha = "alpha_" + str(number)
        gamma = "gamma_" + str(number)
        if alpha in alpha_genes or gamma in gamma_genes:
            list_genes.append("genes_" + str(number))
    comparisons_genes = list(itertools.combinations(list_genes, 2)) # All possible combinations of epitopes counts in amplicons for pair-wise comparisons.
    dict_criteria2_genes = comparison(comparisons_genes, Epsilon, "genes")
    ahpy_criteria2_genes = ahpy.Compare(name="Genes", comparisons=dict_criteria2_genes, precision=3)
    return(ahpy_criteria2_genes, list_genes)

def criteria3epitope(rawdata, Epsilon, epitopes_with_score):
    genotypes = list(rawdata.iloc[:,0])
    genotypes_pairs = list(itertools.combinations(genotypes, 2))
    new_data = rawdata # Copy.
    new_data = new_data.set_index('genotype')
    list_columns = ["genotype"] + epitopes_with_score
    data_epitopes = rawdata.loc[:,list_columns]
    ahpy_criteria3_epitopes = AHP_preparation(data_epitopes, new_data, genotypes_pairs, Epsilon)
    return(ahpy_criteria3_epitopes)

def criteria3gene(rawdata, Epsilon, list_genes):
    data_alpha_genes = rawdata.iloc[:,np.r_[0,103:110]]
    data_gamma_genes = rawdata.iloc[:,np.r_[0,111:123]]
    data_genes = pd.DataFrame()
    data_genes["genotype"] = data_alpha_genes["genotype"]
    genotypes = list(rawdata.iloc[:,0])
    genotypes_pairs = list(itertools.combinations(genotypes, 2))
    new_data = rawdata # Copy.
    new_data = new_data.set_index('genotype')
    for gene in list_genes:
        number = gene.split("_")[1].strip()
        if "alpha_" + str(number) in data_alpha_genes.columns and "gamma_" + str(number) in data_gamma_genes.columns:
            data_genes["gene_" + str(number)] = data_alpha_genes["alpha_" + str(number)] + data_gamma_genes["gamma_" + str(number)]
            new_data["gene_" + str(number)] = new_data["alpha_" + str(number)] + new_data["gamma_" + str(number)]
        if "alpha_" + str(number) in data_alpha_genes.columns and "gamma_" + str(number) not in data_gamma_genes.columns:
            data_genes["gene_" + str(number)] = data_alpha_genes["alpha_" + str(number)]
            new_data["gene_" + str(number)] = new_data["alpha_" + str(number)]
        if "alpha_" + str(number) not in data_alpha_genes.columns and "gamma_" + str(number) in data_gamma_genes.columns:
            data_genes["gene_" + str(number)] = data_gamma_genes["gamma_" + str(number)]
            new_data["gene_" + str(number)] = new_data["gamma_" + str(number)]
    ahpy_criteria3_genes = AHP_preparation(data_genes, new_data, genotypes_pairs, Epsilon)
    return(ahpy_criteria3_genes)