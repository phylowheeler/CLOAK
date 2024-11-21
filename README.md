# CLOAK
Repository for the multiple sequence alignment filtering program: Cleaning on Alignment (K)onsensus (CLOAK)

## Usage
This tool can be run by simply downloading the cloak.py file from this repository, and running it with python3. The user must also specify an Ensemble FASTA (EFA) file containing a set of amino acid multiple sequence alignments to be used as the input alignment set for this tool. The resutlting filtered multiple sequence alignment will be output as myfile.cloak.fa
```
python3 cloak.py -alignments myfile.efa
```
This tool can work with any set of input multiple sequence alignemnts, but it is recommended that the user generates a set of variant alignments using the stratified ensemble option in [muscle5](https://www.drive5.com/muscle/).
```
muscle -align sequences.fasta -stratified -output ensemble.efa
```

## Directories in this Repository

### Benchmarking
Scripts to score the performance of multiple sequence alignment filtering programs

### Substitution_Models
Guide for training amino acid substitution models on filtered multiple sequence alignments, and scripts to compare substituion models to each other.

### Figures
Scripts to generate publication figures

## Dependencies

### Python3
This program requires Python3 to run. Python3 can be installed with [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or your prefered method.

### R
There are scripts in this repository that make use of [R](https://www.r-project.org/)
