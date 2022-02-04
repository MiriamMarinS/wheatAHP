import ahpy
import numpy as np
import pandas as pd
from pandas import read_csv
import argparse
from Source.epitope import Epitope
from Source.AHPmethod import criteria1, criteria2epitope, criteria2gene, criteria3epitope, criteria3gene

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "/path/to/input file with usearch frquencies results.")
parser.add_argument("-o", "--path_output", help = "/path/to/ output files directory.", default="./")
parser.add_argument("-e", "--epitopes", help = "/path/to/epitope fasta file.")
parser.add_argument("-d", "--remove_epitope_duplication", help = "Remove epitope with duplicated sequences. Optios = yes | no.", default="yes")
parser.add_argument("-ol", "--oligopeptides", help = "/path/to/oligopeptides tsv file.")
args = parser.parse_args()

def main():
    rawdata = read_csv(args.input, sep = ";")
    Epsilon = 1 # Value to add to frequencies in NGS data to perform Log2(x + Epsilon) for avoid errors on cases when x = 0.
    ahpy_dict = {} # Dict with results for aphy method changing the intensity value between epitopes and genes in criteria1.
    for intensity in range(1, 9 + 1):
        # Criteria 1:
        epitopes_genes_comparisons = {('Epitopes', 'Genes'): intensity}
        ahpy_criteria1 = criteria1(epitopes_genes_comparisons)
        
        #Criteria 2:
        epitopes = Epitope(args.epitopes, args.remove_epitope_duplication) # Epitope sequences library.
        oligoEpitopes = epitopes.findEpoligos(args.oligopeptides) # Find epitopes in 20-amino acid oligopeptides and keep the max score for each one.
        ahpy_criteria2_epitopes, epitopes_with_score = criteria2epitope(oligoEpitopes, Epsilon) # Intensity values for each epitope comparison based on score from oligopeptides.
        ahpy_criteria2_genes, list_genes = criteria2gene(rawdata, Epsilon) # Intensity values for each amplicon comparison based on number of epitopes in them.

        # Criteria 3:
        ahpy_criteria3_epitopes = criteria3epitope(rawdata, Epsilon, epitopes_with_score) # Intensity values between lines for each epitope based on the epitope frequency in each one.
        ahpy_criteria2_epitopes.add_children(ahpy_criteria3_epitopes)
        ahpy_criteria3_genes = criteria3gene(rawdata, Epsilon, list_genes) # Intensity values between lines for each amplicon type based on the amplicon frequency in each one.
        ahpy_criteria2_genes.add_children(ahpy_criteria3_genes)
        ahpy_criteria1.add_children([ahpy_criteria2_epitopes, ahpy_criteria2_genes])
        ahpy_dict[intensity] = ahpy_criteria1.target_weights # Scores for each line.
    
    scores = pd.DataFrame(data=ahpy_dict)
    scoresdata = pd.concat([scores, scores.mean(axis=1).rename('mean'), scores.std(axis=1).rename('sd')], axis=1)
    scoresdata.to_csv(args.path_output + 'table_scores.txt', header=True, sep='\t', mode='a')

if __name__ == "__main__":
    main()
