#!/usr/bin/env python3

import os
from pathlib import Path
from Bio import SeqIO

def fasta_to_nexus(input_dir, output_dir):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for fasta_file in input_dir.glob("*.fa*"):  # matches .fa, .fasta, etc.
        records = list(SeqIO.parse(fasta_file, "fasta"))
        if not records:
            continue

        for r in records:
            r.annotations["molecule_type"] = "protein"

        out_file = output_dir / (fasta_file.stem + ".nex")
        SeqIO.write(records, out_file, "nexus")
        print(f"{fasta_file.name} -> {out_file.name}")

if __name__ == "__main__":
    # Change these paths as needed
    input_folder = "alignments_renamed"
    output_folder = "alignments_nexus"
    fasta_to_nexus(input_folder, output_folder)
