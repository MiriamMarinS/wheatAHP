import math
import ahpy

def comparison(comparisons_list, Epsilon, type, oligoEpitopes=None):
    intensity_high = [9, 8, 7, 6, 5, 4, 3, 2, 1]
    intensity_low = [1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9]   
    higher_ratio = 0.0
    lower_ratio = 0.0
    dict_comparisons = {}
    for comparison in comparisons_list:
        if type == "epitopes":
            x = oligoEpitopes[comparison[0]]
            y = oligoEpitopes[comparison[1]]
        elif type == "genes":
            x = float(comparison[0].split("_")[1].strip())
            y = float(comparison[1].split("_")[1].strip())
        try:
            #FC: log2(FC) = log2(B) - log2(A), B > A, but in AHP a FC > 0 (> intensity 1) indicates that A > B,
            # so the order is changed in the comparison to: log2(FC) = log2(A) - log2(B).
            z = math.log(x + Epsilon, 2) - math.log(y + Epsilon, 2)
        except:
            continue
        dict_comparisons[comparison] = z
        if z > higher_ratio:
            if z > 0.0:
                higher_ratio = z
        if z < lower_ratio:
            if z < 0.0:
                lower_ratio = z    
    if higher_ratio < 0.0:
        print("The higher ratio is not upper than 1.0.")
    if lower_ratio > 0.0:
        print("The lower ratio is upper than 1.0.")
    new_dict_comparisons = {}

    for k, v in dict_comparisons.items():
        if v > 0.0:
            w = (v*9)/higher_ratio
            new_w = min(intensity_high, key=lambda i:abs(i-w))
        elif v < 0.0:
            #w = (z*(1/9))/lower_ratio
            w = (v*9)/lower_ratio
            w = 1/w
            new_w = min(intensity_low, key=lambda i:abs(i-w))
        else:
            new_w = 1
        new_dict_comparisons[k] = new_w
    return(new_dict_comparisons)

class Linefreq():
    def __init__(self, name, comparisons_dict):
        self.name = name
        self.comparisons_dict = comparisons_dict

def AHP_preparation(data_epitopes, new_data, genotypes_pairs, Epsilon):
    epitopes = list(data_epitopes.columns[1:])
    intensity_high = [9, 8, 7, 6, 5, 4, 3, 2, 1]
    intensity_low = [1, 1/2, 1/3, 1/4, 1/5, 1/6, 1/7, 1/8, 1/9]
    comparisons = []
    dict_epitopes = {}
    for epitope in epitopes:
        comparison_dict = {}
        higher_ratio = 0.0
        lower_ratio = 0.0
        for genotypes in genotypes_pairs:
            try:
                ratio = math.log(float(new_data[epitope][genotypes[0]]) + Epsilon, 2) - math.log(float(new_data[epitope][genotypes[1]]) + Epsilon, 2)
            except:
                continue
            comparison_dict[genotypes] = ratio
            if ratio > higher_ratio:
                if ratio > 0.0:
                    higher_ratio = ratio
            if ratio < lower_ratio:
                if ratio < 0.0:
                    lower_ratio = ratio
        if higher_ratio < 0.0:
            print("The higher ratio is not upper than 0.")
        if lower_ratio > 0.0:
            print("The lower ratio is upper than 0.")
        new_comparison_dict = {}
        for k, v in comparison_dict.items():
            if v == 0.0:
                new_comparison_dict[k] = 1
            if v > 0.0:
                value = (v*9)/higher_ratio
                new_value = min(intensity_high, key=lambda x:abs(x-value))
                new_comparison_dict[k] = new_value
            if v < 0.0:
                value = (v*9)/lower_ratio
                value = 1/value
                new_value = min(intensity_low, key=lambda x:abs(x-value))
                new_comparison_dict[k] = new_value
        dict_epitopes[epitope] = Linefreq(epitope, new_comparison_dict)
        comparisons.append(ahpy.Compare(epitope, new_comparison_dict, precision=3))
    return(comparisons)