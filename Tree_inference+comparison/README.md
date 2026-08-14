Scripts used to generate phylogenetic trees, and to perform calculations using the resulting trees.
- Astral.sh: used to calculate quartet support. User must specify file directory of gene trees and file path to species tree
- calculate_lrm_distance.py: calculates Lin-Rajan-Moret distance between a gene tree and a species tree. Specify both trees as parameters
- fasta_to_nexus.py: converts fasta alignment to nexus alignment
- mcoffee.sh: used to find a consensus tree from a set of input trees. User must specify the path to a directory containing the input trees
- run_lrm_distance.sh: calculates the Lin-Rajan-Moret distance between all trees in a directory and a species tree, both given as parameters
- tree_inference.sh: script to run IQ-Tree. User must specify an input directory of alignments, desired output directory, and the empirical substitution model desired. The script will automatically use empirical amino acid frequencies (+f) and choose a rate-heterogeneity among sites model using ModelFinder. To run IQ-Tree with a GTRpMIX profile mixture model, add "+C60" to the end of the substitution model name, and add the "--link-exchange" flag to the run. 
