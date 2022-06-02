from genericpath import exists
from Bio import SeqIO
import pandas as pd
import numpy as np
from Source.utils import variants
import math

class Epitope():
    def __init__(self, epitopesFile, removeDuplications):
        self.epitopesFile = epitopesFile
        self.removeDuplications = removeDuplications
        self.library = self.epitope_library()
        if self.removeDuplications == "yes":
            self.library = self.epitope_library_remove_duplication()

    def epitope_library(self):
        self.epitope_sequences = {}
        for record in SeqIO.parse(self.epitopesFile, "fasta"):
            self.epitope_sequences[str(record.id)] = str(record.seq)
        return(self.epitope_sequences)

    def epitope_library_remove_duplication(self):
        # In the case of original, deaminated1, deaminated2 and moAb epitopes methods, the epitopes are grouped in classes.
        # If epitopes in the same class have the same sequences, they will be considered as one epitope. In other methods, the duplication will be removed in the group of all epitopes, because they are not grouped in classes.

        #Find duplicated sequences in epitope library.
        duplicated_sequences = {}
        for epitope_id, value in self.library.items():
            epitope_seq = value.strip()
            duplicated_sequences.setdefault(epitope_seq, []).append(epitope_id)

        #New dictionary with epitopes grouped by their sequences separated in original, deaminated1 and deaminated2 groups, or unclassified group for the other methods.
        self.new_epitope_library = {}
        for k, v in duplicated_sequences.items():
            self.new_epitope_library[';'.join(v)] = k

        #The name of epitopes with duplicated sequences are concatenated by ;. But the name is too long. So we kept the first epitope name, and it can be checked in epitope library printed.
        self.epitope_sequences2 = {}
        for k, v in self.new_epitope_library.items():
            if len(k.split(";")) > 1:
                self.epitope_sequences2[k.split(";")[0]] = v
            else:
                self.epitope_sequences2[k] = v
        return(self.epitope_sequences2)

    def findEpoligos(self, oligosFile):
        self.oligosFile = oligosFile
        peptidesData = pd.read_csv(self.oligosFile, sep = "\t", header = 0)
        epitopes_find = {}
        for epitopeID, epitopeSeq in self.library.items():
            for index, row in peptidesData.iterrows():    
                count = 0 # Takes values 0 | 1, if epitope in oligopeptide.
                if epitopeSeq in row["20-amino acid oligopeptide"]: count = 1
                epitopes_find.setdefault(epitopeID, []).append(count)
        peptidesEpitopes = pd.DataFrame.from_dict(epitopes_find)
        new_peptidesData = pd.merge(peptidesData, peptidesEpitopes, left_index=True, right_index=True)

        # Keep the max score that an epitope takes in all the peptides.
        epitopeScores = {}
        for index, row in new_peptidesData.iterrows():
            for epitopeID, epitopeSeq in self.library.items():
                if row[epitopeID] > 0:
                    epitopeScores.setdefault(epitopeID, []).append(row["First-round score"])
        epitopemaxScore = {}
        for epitopeID, scores in epitopeScores.items():
            epitopemaxScore[epitopeID] = max(scores)
        return(epitopemaxScore)
    
    def epitopemismatch(self, mismatch_allowed, table_freq, amp_dict):
        self.amp_dict = amp_dict
        self.mismatch_allowed = mismatch_allowed
        self.table_freq = table_freq
        index = 0
        for column in list(table_freq.columns): # column = line.
            if index == 0:
                self.table_temp = variants(self.library, self.amp_dict, self.mismatch_allowed, self.table_freq, column)
            else:
                self.table_temp = pd.merge(self.table_temp,
                    variants(self.library, self.amp_dict, self.mismatch_allowed, self.table_freq, column)
                    , left_on = ['Variant', 'Original'], right_on = ['Variant', 'Original'], how = 'outer')
            index += 1
        return(self.table_temp)
    
    def epitope_index(self, table_epitope_freq):
        self.table_epitope_freq = table_epitope_freq
        dict_epitopes = {}
        epsilon = 0.0000001
        dict_epitopes["Variant"] = self.table_temp["Variant"]
        dict_epitopes["Original"] = self.table_temp["Original"]
        for index, row in self.table_temp.iterrows():
            for column in list(self.table_temp.columns)[2:]:
                if row["Original"] in self.table_epitope_freq.columns:
                    freq_original = self.table_epitope_freq[row["Original"]][self.table_epitope_freq.index.values == column].values[0]
                    freq_variant = row[column]
                    ep = math.log(freq_variant + epsilon, 2) - math.log(freq_original + epsilon, 2)
                    #ep = freq_variant/(freq_original + epsilon)
                    #ep = row[column]/(self.table_epitope_freq[row["Original"]][self.table_epitope_freq.index.values == column].values[0] + epsilon)
                    dict_epitopes.setdefault(column, []).append(ep)
                else:
                    dict_epitopes.setdefault(column, []).append(np.NaN)
        self.table_epitope = pd.DataFrame(dict_epitopes, columns = list(self.table_temp.columns)) #table_variants = pd.DataFrame([k + (v,) for k, v in dict_variants.items()], columns=['Variant', 'Original', column]))
        self.table_epitope = self.table_epitope.dropna()
        self.table_epitope = self.table_epitope.reset_index(drop=True)
        return(self.table_epitope)
