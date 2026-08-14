#!/bin/bash

module load anaconda
source ~/.bashrc && conda activate
conda activate t_coffee

for file in /path/to/directory/*; 
do [ -f "$file" ] && t_coffee "$file" -mode mcoffee ; done