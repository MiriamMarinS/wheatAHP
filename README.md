# wheatAHP
Implementation of the Analytic Hierarchy Process (AHP) to obtain the immunogenicity score of wheat lines based on their celiac disease (CD) epitopes matches on alpha- and gamma-gliadins amplicons by NGS and the score for oligopeptides based on an IFN-g ELISpot assays with fresh peripherical blood mononuclear cells (PBMCs).

# **Immunogenic score for 20-amino acid oligopeptides of wheat**

In the work of Tye-Din *et al.* (2010), they present a table of content with a score for each 20–amino acid oligopeptide based on their average relative frequency (mean normalized response of donor who responded to at least one peptide) of specific T cells present in blood. The 20-amino acid oligopeptides were pretreated with tissue transglutaminase (tTG). The wheat/barley/rye gluten peptide library design is described in Tye-Din's work.
\
\
\
*TYE-DIN, Jason A., et al. Comprehensive, quantitative mapping of T cell epitopes in gluten in celiac disease. Science translational medicine, 2010, vol. 2, no 41, p. 41ra51-41ra51.*

# **Search epitopes in oligopeptides**
We search for CD epitopes in the peptide library.
The 20-amino acid oligopeptides related to scores were presented without deamidation: we search non-deamidated epitopes in the peptide library. We assign to each epitope the maximum score value of the oligopeptide in which it is found.

# **AHP method**
With the ahpy python module (```pip install ahpy```), we assign a immunogenicity score to each genotype.
The method is divided on three criteria:
* Criteria 1: Epitopes vs amplicons.
* Criteria 2:
  * Criteria 2 epitopes:
    * Pair-wise comparisons between epitopes based on Tye-Din immunogenicity scores.
  * Criteria 2 amplicons:
    * Pair-wise comparisons between amplicons based on the number of epitopes in each one.
* Criteria 3:
  * Criteria 3 epitopes.
    * Pair-wise comparisons between lines based on the frequency of each epitope in them.
  * Criteria 3 amplicons.
    * Pair-wise comparisons between lines based on the frequency of each amplicon type (by number of epitopes in them) in each one.

![alt text](./results/Diagram.png?raw=true)

```mermaid
flowchart LR;

  A(Diagram) ---> B(Amplicons);
  A ---> C(Epitopes);
  B ---> A0(Amp_0) & A1(Amp_1) & A2(Amp_2) & A3(Amp_3) & A4(Amp_4) & A5(Amp_5) & A6(Amp_6) & A7(Amp_7) & A8(Amp_8) & A9(Amp_9) & A10(Amp_10) & A11(Amp_11) & A12(Amp_12) & A13(Amp_13) & A14(Amp_14);
  C ---> E0(DQ2.5_glia_a1a) & E1(DQ2.5_glia_a1b) & E2(DQ2.5_glia_a2) & E3(DQ2.5_glia_g2a) & E4(DQ2.5_glia_g4e) & E5(DQ2.5_glia_g4d) & E6(DQ2.5_glia_g4b) & E7(DQ2.5_glia_g5) & E8(DQ2.5_glia_g4c) & E9(DQ2.5_glia_g3) & E10(DQ2.5_glia_g4a) & E11(DQ2.5_glia_w1) & E12(DQ2.5_glia_hor_2);
```

**Figure.** Diagram of the three criteria for AHP method.

# **Running the method**

```
python wheatAHP.py -i ./Input/Data_matrix_AHP.csv -o ./results/ -e ./Input/epitopes.fasta -d yes -ol ./Input/20_oligopeptides_scores.txt
```

# **Scores and AHP graph**
For each genotypes, a immunogenic score is assigned. For representation, the mean of score values changing the instensity value between epitopes and amplicons pair-wise comparisons in criteria 1 is calculated, and the standar deviation.

*Output: table_scores.txt*

Figure_scores.R

![alt text](./results/Figure_scores.png?raw=true)

**Figure.** Immunogenicity score per each genotype. The genus (BW: bread wheat, DW: durum wheat, HT: Tritordeum) and the presence of rye (N: no, Y: yes) are indicated. The mean of scores (changing intensity values in the criteria 1) and the standard deviation (bars) are represented.

# **Searching epitopes in NGS amplicons**
For searching epitopes and epitope variants with 1 or more mismatches in amplicon sequences of alpha and gamma gliadins and calculating their relative abundance in lines. In addition, calculation of the ratio of relative abundance variation between the epitope variant and the original epitope for each line is also implemented.

```
python ./findEpmismatch.py -i </path/to/file with amplicon frequencies in lines> -o </path/to/output dir> -e </path/to/epitopes fasta file> -d <remove epitope duplications: yes | no> -p <output prefix> -m <number of mismatches: 0, 1, ...> -a </path/to/amplicon peptides fasta file> -ep </path/to/file with epitope frequencies in lines>
```

**Output:**
* File with frequencies of epitope variants in lines, tsv format.
* File with ratio of frequency variation between the epitope variant and the original epitope, tsv format.