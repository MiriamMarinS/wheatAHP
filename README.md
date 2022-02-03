# wheatAHP
Implementation of the Analytic Hierarchy Process (AHP) to obtain the immunogenicity score of wheat lines based on their celiac disease (CD) epitopes matches on alpha- and gamma-gliadins amplicons by NGS and the score for oligopeptides based on an IFN-g ELISpot assays with fresh peripherical blood mononuclear cells (PBMCs).

# **Immunogenic score for 20-amino acid oligopeptides of wheat**

In the work of Tye-Din *et al.* (2010), they present a table of content with a score for each 20–amino acid oligopeptide based on their average relative frequency (mean normalized response of donor who responded to at least one peptide) of specific T cells present in blood. The 20-amino acid oligopeptides were pretreated with tissue transglutaminase (tTG). The wheat/barley/rye gluten peptide library design is described in Tye-Din's work.

We search for CD epitopes in the peptide library.
The 20-amino acid oligopeptides related to scores were presented without deamidation: we search non-deamidated epitopes in the peptide library.

TYE-DIN, Jason A., et al. Comprehensive, quantitative mapping of T cell epitopes in gluten in celiac disease. Science translational medicine, 2010, vol. 2, no 41, p. 41ra51-41ra51.

```
python wheatAHP.py -i ./Input/Data_matrix_AHP.csv -o ./results/ -e ./Input/epitopes.fasta -d yes -ol ./Input/20_oligopeptides_scores.txt
```

# **Scores**
For each genotypes, a immunogenic score is assigned.
*Output: table_scores.txt*
![alt text](https://github.com/MiriamMarinS/wheatAHP/results//Figure_score.png?raw=true)