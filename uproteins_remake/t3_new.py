"""gather new t3 sequences that werent in the previous t3 version of analyses of smorfs for uproteins paper on mtb"""
from Bio import SeqIO

old_smorfs = "/home/adriana/tcc/smorfs/sequences/gorfs_T4.fasta"
new_smorfs = "/home/adriana/uproteins/redo/240515_smorfTiers_genome.txt"
output_file = "/home/adriana/uproteins/redo/new_T3_gorfs_2.txt"

# Read and store all old smorfs
old_smorfs_set = set(record.description for record in SeqIO.parse(old_smorfs, 'fasta'))

# Initialize an empty list to store new T3 smorfs
new_smorfs_list = []

# Read new smorfs file and find new T3 smorfs that are not in the old smorfs set
with open(new_smorfs, 'r') as smorf_file:
    next(smorf_file)  # Skip the header

    for line in smorf_file:
        smorf, tier = line.strip().split('\t')

        if smorf not in old_smorfs_set and tier == 'T3':
            print(f"new smorf found {smorf}")
            new_smorfs_list.append(smorf)

# Print the new T3 smorfs
print(new_smorfs_list)
#
with open(output_file, 'w') as f:
    for smorf in new_smorfs_list:
        f.write(f"{smorf}\n")

print(f"New T3 smorfs have been saved to {output_file}")







