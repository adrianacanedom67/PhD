"""gather new T3 sequences from csv file"""

import pandas as pd

# df = pd.read_csv('/home/adriana/uproteins/redo/240516_cat_genomeTranscriptome_pepFDR_tiers_filtered_proteins.csv')
csv_path = '/home/adriana/uproteins/redo/240516_cat_genomeTranscriptome_pepFDR_tiers_filtered_proteins.csv'
new_smorfs = '/home/adriana/uproteins/redo/new_T3_torfs.txt'

try:
    df = pd.read_csv(csv_path, on_bad_lines='skip', sep="\t")

except pd.errors.ParserError as e:
    print(f"Error reading the CSV file: {e}")
    # Handle the error or investigate the problematic lines
    df = pd.read_csv(csv_path, on_bad_lines='skip', delimiter=',')

# Read the new smorf names into a list
with open(new_smorfs, 'r') as smorf_file:
    new_smorfs_list = [line.strip() for line in smorf_file if line.strip()]

# print(new_smorfs_list)


smorf_dict = {} # store new smorfs (keys) and their corresponding sequences (values)

for index, row in df.iterrows():
    smorf_names = row['Final Entries']
    sequences = row['Extended Sequence']
    if smorf_names in new_smorfs_list:
        smorf_dict[smorf_names] = sequences
print(smorf_dict)

fasta_output_path = "/home/adriana/uproteins/redo/new_T3_torfs.fasta"

with open(fasta_output_path, 'w') as fasta_file:
    for smorf_name, sequence in smorf_dict.items():
        fasta_file.write(f'>{smorf_name}\n')
        fasta_file.write(f'{sequence}\n')

print(f"Fasta file successfully generated and saved to {fasta_output_path}")



