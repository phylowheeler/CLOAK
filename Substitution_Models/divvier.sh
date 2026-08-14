#!/bin/bash

module load anaconda
source ~/.bashrc && conda activate
conda activate divvier

for file in /path/to/directory/*; 
do [ -f "$file" ] && divvier -divvygap "$file"  ; done