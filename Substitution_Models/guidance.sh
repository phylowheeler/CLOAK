#!/bin/bash

for file in /path/to/directory/*; 
do [ -f "$file" ] && 

perl guidance.pl --seqFile "$file" --msaProgram MUSCLE --seqType aa --outDir {output directory}; done
