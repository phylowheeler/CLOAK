#!/bin/bash

for file in /path/to/directory/*; 
do [ -f "$file" ] && 
julia correction_multi.jl  -m - "$file" ; done