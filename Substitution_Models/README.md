Directory containing methods to generate and evaluate amino acid substitution models trained on filtered multiple sequence alignments.

## Model Training
Amino Acid substitution models were trained using Qmaker, implemented in [iqtree2](http://www.iqtree.org/) as described [here](http://www.iqtree.org/doc/Estimating-amino-acid-substitution-models#estimating-a-model-from-a-single-concatenated-alignment). Tree inference was also performed using iqtree2.

## Models
Directory containing the substitution model files

## Trees
Directory containing the raw tree files

### matrices.R
Script to read in and compare substitution models

### mutational_accessibility.csv
table containing the data on exchangeability parameters according to mutational accessibility

### treedistance.py
Script for calculating the lin-rajan-moret distance between phylogenetic trees

### tree_distance.csv
table containing tree distances using different substitution models
