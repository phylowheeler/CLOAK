# CLOAK
Repository for the multiple sequence alignment filtering program: Cleaning on Alignment (K)onsensus (CLOAK). This software tool is designed to filter out errors from amino acid multiple sequence alignments by identifying dissimmilarities between variant alignments.

## Usage
This tool can be used by simply downloading the cloak.py file from this repository, and running it with python3. The user must also specify an Ensemble FASTA (EFA) file containing a set of amino acid multiple sequence alignments to be used as the input alignment set for this tool. The resutlting filtered multiple sequence alignment will be output as myfile.cloak.fa
```
python3 cloak.py -alignments myfile.efa
```
This tool can work with any set of input multiple sequence alignemnts, but it is recommended that the user generates a set of variant alignments using the stratified ensemble option in [muscle5](https://www.drive5.com/muscle/).
```
muscle -align sequences.fasta -stratified -output ensemble.efa
```
CLOAK can also be used directly within Muscle5 with the -cloak option, using providing either a multiple alignment file or directory of variant alignments.
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
