import glob
from Bio import SeqIO, AlignIO
import os
import csv

def get_fasta_sequence_lengths(directory_path):
	fasta_files = glob.glob(os.path.join(directory_path, '*.fasta')) + \
					glob.glob(os.path.join(directory_path, '*.fa')) + \
					glob.glob(os.path.join(directory_path, '*.fna')) + \
					glob.glob(os.path.join(directory_path, '*.aln')) 
	if not fasta_files:
		print(f"No FASTA files found in the directory: {directory_path}")
		return
	print(f"Found {len(fasta_files)} FASTA files.")
	lengths = []
	for file_path in fasta_files:
		file_name = os.path.basename(file_path)
		try:
			alignment = SeqIO.parse(file_name, "fasta")
			align_length = []
			for record in alignment:
				rec_length = len(str(record.seq).replace("-", ""))
				align_length.append(rec_length)
			mean_length = sum(align_length) / len(align_length)
			lengths.append(mean_length)
		except Exception as e:
			print(f"Error processing file {file_name}: {e}")
	return lengths


lengths=get_fasta_sequence_lengths(".")

with open('mean_lengths.csv', 'w', newline='') as f:
	writer = csv.writer(f)
	for item in lengths:
		writer.writerow([item])