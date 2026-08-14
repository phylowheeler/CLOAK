Directory containing methods to generate and evaluate amino acid substitution models trained on filtered multiple sequence alignments.

## Alignment filtering
Prior to substitution model training, all datasets were realigned using Muscle5. Filtering was then performed using a variety tools as described in the paper methods

## Model Training
Amino Acid substitution models were trained using Qmaker, implemented in [iqtree2](http://www.iqtree.org/) as described [here](http://www.iqtree.org/doc/Estimating-amino-acid-substitution-models#estimating-a-model-from-a-single-concatenated-alignment). Tree inference was also performed using iqtree2. Instructions can also be seen in the qmaker.sh file

## Pretrained_Q
Directory containing the substitution model files
