"""generate fasta files for all new tiers up to T3"""

from Bio import SeqIO
import pandas as pd

tier_information = '/home/adriana/uproteins/redo/240515_smorfTiers_transcriptome.txt'
sequences_path = '/home/adriana/uproteins/redo/240516_cat_genomeTranscriptome_pepFDR_tiers_filtered_proteins.csv'
old_smorfs = '/home/adriana/tcc/smorfs/sequences/torfs_T1.fasta'

new_smorf_list = [] # new T1_smorfs from tier_information
current_smorf_list = [] # old T1_smorfs from old_smorfs that are also in current new tier
old_smorf_list = [] # all the old T1_smorfs
discarded_smorf_list = [] # old T1_smorfs from old_smorfs that arent in tier_information anymore and will be discarded

records = SeqIO.parse(old_smorfs, 'fasta')
for record in records:
    old_smorf = record.description
    old_smorf_list.append(old_smorf)

# print(old_smorf_list)

with open(tier_information, 'r') as tier_file:
    next(tier_file)
    for line in tier_file:
        smorf, tier = line.strip().split('\t')

        if smorf not in old_smorf_list and tier == 'T1':
            new_smorf_list.append(smorf)
            current_smorf_list.append(smorf)
            print(f'new smorf {smorf} found')
        elif smorf in old_smorf_list and tier =='T1':
            print(f'smorf {smorf} already in old Fasta file')
            current_smorf_list.append(smorf)

# Check which old smorfs are no longer in the new tier information
for old_smorf in old_smorf_list:
    if old_smorf not in current_smorf_list:
        print(f'Smorf {old_smorf} not in new Tier classification')
        discarded_smorf_list.append(old_smorf)



print(new_smorf_list)
old_smorf_list = list(set(old_smorf_list))
print(discarded_smorf_list)

smorf_dict = {}  # dictionary that will store the current smorf names as keys
# and their corresponding sequences as values

df = pd.read_csv(sequences_path, on_bad_lines='skip', sep='\t')
for index, row in df.iterrows():
    smorf_names = row['Final Entries']
    sequences = row['Extended Sequence']
    if smorf_names in current_smorf_list:
        smorf_dict[smorf_names] = sequences

print(smorf_dict)

fasta_output_path = "/home/adriana/uproteins/redo/current_T1_torfs.fasta"

with open(fasta_output_path, 'w') as fasta_file:
    for smorf_name, sequence in smorf_dict.items():
        fasta_file.write(f'>{smorf_name}\n')
        fasta_file.write(f'{sequence}\n')

print(f"Fasta file successfully generated and saved to {fasta_output_path}")
