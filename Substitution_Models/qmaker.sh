#!/bin/bash


# Check for input argument
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch Q_maker.sh <model_name>"
    exit 1
fi

model_name=$1

module load anaconda
source ~/.bashrc && conda activate
conda activate iqtree3

#From a folder
# step 1: infer a  separate tree for each alignment with reversible models as initial models
# iqtree3 -seed 1 -T AUTO -safe -S surface_alignments_nexus -mset LG,WAG,JTT -cmax 4 -pre surface_alignments_nexus

# # step 2: estimate a join reversible matrix across all alignments
# iqtree3 -seed 1 -T AUTO -S surface_alignments_nexus.best_model.nex -te surface_alignments_nexus.treefile --model-joint GTR20+FO --init-model LG -pre surface_alignments_nexus.GTR20

# # step 3: extract the resulting reversible matrix
# grep -A 21 "can be used as input for IQ-TREE" surface_alignments_nexus.GTR20.iqtree | tail -n20 > Q.surface_alignments_nexus

iqtree3 --seed 1 -T AUTO -s ${model_name}.nex -m MFP -mset LG,WAG,JTT -cmax 4 --prefix ${model_name}

iqtree3 -seed 1 -T AUTO -s ${model_name}.nex -te ${model_name}.treefile --init-model LG --model-joint GTR20+FO --prefix ${model_name}_GTR20_FO

grep -A 21 "can be used as input for IQ-TREE" ${model_name}_GTR20_FO.iqtree | tail -n20 > Q.${model_name}
