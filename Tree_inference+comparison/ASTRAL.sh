#!/bin/bash

module load anaconda
source ~/.bashrc && conda activate


for file in /path/to/directory/*; 
do [ -f "$file" ] && 
java -jar astral.5.7.1.jar -i "$file" -q {species treefile}
; done
