# CLOAK
Repository for the multiple sequence alignment filtering program: Cleaning on Alignment (K)onsensus (CLOAK). This software tool is designed to filter out errors from amino acid multiple sequence alignments by identifying dissimmilarities between variant alignments.

## Usage
This tool can be used in one of two ways. First, the python version, cloak.py, is available for download from this repository. Alternatively, a version of CLOAK has been integrated directly into [muscle5] (https://www.drive5.com/muscle/). Instructions for running both versions are provided below. We recommend using the muscle implementation, as it has additional functionality and is will be better supported.

For both versions, the user must provide a set of multiple sequence alignments as input, either as an Ensemble FASTA (EFA) file, or the path to a directory containing the multiple sequence alignment files in FASTA format. This tool can work with any set of input multiple sequence alignments. It has been tested with the 16 variant alignments inferred with muscle5 using the stratified ensemble option:
```
muscle -align sequences.fasta -stratified -output ensemble.efa
```
The python version can be run in python3:
```
python3 cloak.py -alignments myfile.efa
```
The muscle version can be run as a separate command within muscle:
```
muscle -cloak input_ensemble_file -mincol <integer> -output <output_file_name>
```
Arguments (for muscle version):
- -input_ensemble_file : Path to the input MSA file, which can either be an EFA file 
                        or a text file with paths to individual MSAs on each line
- -mincol <integer>   : Minimum number of non-gap characters required per column
                        for that column to be retained in the output.
                        Default value of 2 if not specified
- -output <filename>  : Name of the file where the filtered MSA will be written. By default this will be {input_file_name}.cloak.fa

## Directories in this Repository

### MSA_benchmarking
Scripts to score the performance of multiple sequence alignment filtering programs based on the BAliBASE benchmark

### Substitution_Models
Guide and scripts for training amino acid substitution models on filtered multiple sequence alignments. Pre-trained substitution models. Scripts may require editing e.g. to set file paths or HPC settings

### Tree_inference+comparison
Scripts to infer phylogenetic trees and compare them to a reference species tree. Scripts may require editing e.g. to set file paths or HPC settings

### Figures
Scripts and datafiles to generate publication figures

### Data
Raw datafiles used in analyses

## Dependencies
The software environment in the cloak_env.yml file contains the dependencies needed to run the software in this repository. Set up the environment with conda using
```
conda env create --file cloak_env.yml
conda activate cloak_env
```
## License

This code is released under the GNU General Public License v3.0. See [`LICENSE`](LICENSE).
