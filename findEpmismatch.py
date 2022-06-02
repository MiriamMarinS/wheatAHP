import pandas as pd
import numpy as np
from Source.epitope import Epitope
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help = "/path/to/input file with usearch frquencies results of amps.")
parser.add_argument("-o", "--path_output", help = "/path/to/ output files directory.", default="./")
parser.add_argument("-e", "--epitopes", help = "/path/to/epitope fasta file.")
parser.add_argument("-d", "--remove_epitope_duplication", help = "Remove epitope with duplicated sequences. Optios = yes | no.", default="yes")
parser.add_argument("-p", "--prefix", help = "Prefix for the output name.")
parser.add_argument("-m", "--mismatch", help = "Mismatch allowed in epitope search: int.")
parser.add_argument("-a", "--amp_fasta", help = "/path/to/file with amp seq.")
parser.add_argument("-ep", "--epitope_freq", help = "/path/to/file with epitope freq in lines.")
args = parser.parse_args()

table_freq = pd.read_csv(args.input, sep = "\t", header = 0, index_col = 0) # Table of freq of amps.
amp_dict = {str(record.id):str(record.seq) for record in SeqIO.parse(args.amp_fasta, "fasta")} # Dict with amp sequences.
epitopes = Epitope(args.epitopes, args.remove_epitope_duplication) # Epitope sequences library.
data_mismatch = epitopes.epitopemismatch(int(args.mismatch), table_freq, amp_dict)
data_mismatch.to_csv(args.path_output + args.prefix + '_table_mismatch.txt', header=True, sep='\t', mode='a')

table_freq_epitopes = pd.read_csv(args.epitope_freq, sep = "\t", header = 0, index_col = 0)
data_index = epitopes.epitope_index(table_freq_epitopes)
data_index.to_csv(args.path_output + args.prefix + '_table_mismatch_index.txt', header=True, sep='\t', mode='a')