#!/bin/bash


# Check for input argument
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch Q_maker.sh <model_name>"
    exit 1
fi

model_name=$1

module load anaconda
source ~/.bashrc && conda activate
conda activate iqtree2

# step 1: infer a  separate tree for each alignment with reversible models as initial models
iqtree3 -seed 1 -T AUTO -S ${model_name}.nex -mset LG,WAG,JTT -cmax 4 -pre ${model_name}.nex

# step 2: estimate a join reversible matrix across all alignments
iqtree3 -seed 1 -T AUTO -S ${model_name}.best_model.nex -te ${model_name}.treefile --model-joint GTR20+FO --init-model LG -pre ${model_name}.GTR20

# step 3: extract the resulting reversible matrix
grep -A 21 "can be used as input for IQ-TREE" ${model_name}.GTR20.iqtree | tail -n20 > Q.${model_name}
