Directory containing methods to generate and evaluate amino acid substitution models trained on filtered multiple sequence alignments.

## Alignment filtering
Prior to substitution model training, all datasets should be realigned using Muscle5. The muscle alignments should then be filtered with the desired method. For all alignment filtering scripts, users must specify the directory of the input alignments. Guidance2 is set to runs Muscle5 automatically, so it should be given unaligned sequences as inputs.

## Model Training
Amino Acid substitution models were trained using Qmaker, implemented in [iqtree2](http://www.iqtree.org/) as described [here](http://www.iqtree.org/doc/Estimating-amino-acid-substitution-models#estimating-a-model-from-a-single-concatenated-alignment). Instructions can also be seen in the qmaker.sh file

## Pretrained_Q
Directory containing the substitution model files
