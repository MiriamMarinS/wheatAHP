from Bio import SeqIO
import pandas as pd
import numpy as np

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