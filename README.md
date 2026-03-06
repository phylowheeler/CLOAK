# CLOAK
Repository for the multiple sequence alignment filtering program: Cleaning on Alignment (K)onsensus (CLOAK). This software tool is designed to filter out errors from amino acid multiple sequence alignments by identifying dissimmilarities between variant alignments.

## Usage
This tool can be used in one of two ways. First, the python version, cloak.py, is available for download from this repository. Alternatively, a version of CLOAK has been integrated directly into Muscle5. Instructions for running both versions are provided below. For both versions, the user must provide a set of multiple sequence alignments as input, either as an Ensemble FASTA (EFA) file, or the path to a directory containing the multiple sequence alignmet files in FASTA format. 

The python version can be run in python3 as shown below. The resutlting filtered multiple sequence alignment will be output as myfile.cloak.fa
```
python3 cloak.py -alignments myfile.efa
```
This tool can work with any set of input multiple sequence alignemnts, but it is recommended that the user generates a set of variant alignments using the stratified ensemble option in [muscle5](https://www.drive5.com/muscle/), as shown below. 
```
muscle -align sequences.fasta -stratified -output ensemble.efa
```
The muscle version can be run as a separate command within muscle as shown below.
```
muscle -cloak input_ensemble_file -mincol <integer> -output <output_file_name>
```
Arguments:
- input_ensemble_file : Path to the input MSA file, which can either be an EFA file 
                        or a text file with paths to individual MSAs on each line
- -mincol <integer>   : Minimum number of non-gap characters required per column
                        for that column to be retained in the output.
                        Default value of 2 if not specified
- -output <filename>  : Name of the file where the filtered MSA will be written.

## Directories in this Repository

### Benchmarking
Scripts to score the performance of multiple sequence alignment filtering programs

### Substitution_Models
Guide for training amino acid substitution models on filtered multiple sequence alignments, and scripts to compare substituion models to each other.

### Figures
Scripts to generate publication figures

## Dependencies
The software environment in the cloak_env.yml file contains the dependencies needed to run the software in this repository. Set up the environment with conda using
```
conda env create --file cloak_env.yml
conda activate cloak_env
```
