from Bio import SeqIO
import os
import pandas as pd


old_fasta_files_directory = "/home/adriana/uproteins/redo/fasta_files"
fasta_dict = {}
old_fasta_files = os.listdir(old_fasta_files_directory)
# print(old_fasta_files)
for file in old_fasta_files:
    records = SeqIO.parse(f'{old_fasta_files_directory}/{file}', 'fasta')
    # print(f"{file} processed successfuly")
    for record in records:
        identifiers = str(record.description)
        sequence = str(record.seq)
        fasta_dict[identifiers] = [sequence]
        if len(sequence) > 100:
            print(file, identifiers, sequence)

csv_path = "/home/adriana/uproteins/redo/240516_cat_genomeTranscriptome_pepFDR_tiers_filtered_proteins.csv"
df = pd.read_csv(csv_path, on_bad_lines='skip', sep='\t')
smorf_dict = {}
for index, row in df.iterrows():
    smorf_names = row['Final Entries']
    sequences = row['chosen_sequences']
    smorf_dict[smorf_names] = [sequences]


sequences_to_be_changed = []
for key in fasta_dict.keys():
    if key in smorf_dict and fasta_dict[key] != smorf_dict[key]:
        sequences_to_be_changed.append(key)

print(sequences_to_be_changed)


