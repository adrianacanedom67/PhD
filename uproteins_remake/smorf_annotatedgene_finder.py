"""find which novel smorfs overlap with annotated genes"""

import pandas as pd
from Bio import SeqIO

smorfs = '/home/adriana/uproteins/redo/current_T1_gorfs.fasta'
output_file = '/home/adriana/uproteins/redo/T1_gorfs_overlapping_genes.xlsx'

df = pd.read_excel('/home/adriana/uproteins/redo/Supplementary_Table_2-overlapping_genes.xlsx')
pd.set_option('display.max_columns', None)
df = df.fillna(0)
# print(df.head(15))

# this dictionary will store the smorfs (keys) and the corresponding annotated overlapping gene (values)
smorf_annotatedoverlap = {}
all_smorfs = []  # all smorfs in fasta file
overlapped_smorfs = []  # all smorfs in supplementary_table
overlapped_dict = {}   # all smorfs and their overlapping annotated genes in supplementary_table
current_overlappedsmorfs = [] # only smorfs that match in both lists above
current_overlappedsmorfs_dict = {} # dictionary of smorfs that match both with the corresponding genes (values)

for index, row in df.iterrows():
    new_smorf_name = row['novel_ORF']
    overlapping_gene = row['Overlapping gene id'] # "Rv's"
    overlapping_gene_name = row['Overlapping gene name']  # "inha"
    overlapped_smorfs.append(new_smorf_name)
# for new_smorf_name, overlapping_gene in overlapped_smorfs:
    if new_smorf_name not in overlapped_dict:
        overlapped_dict[new_smorf_name] = []  # add more than one value for each key
    overlapped_dict[new_smorf_name].append(overlapping_gene)   # there are more than one overlapping genes for each smorf
    # overlapped_dict[new_smorf_name] = overlapping_gene
# since the original supplementary table has duplicates, remove repeated values inside the same key:
for key in overlapped_dict:
    overlapped_dict[key] = list(set(overlapped_dict[key]))

sorted_dict = {key: overlapped_dict[key] for key in sorted(overlapped_dict)}
# print(sorted_dict)

# for key, value in sorted_dict.items():
#     print(f'{key}: {value}')


# print(f'all in supplementary table: {overlapped_smorfs}')
# overlapped_smorfs = list(set(overlapped_smorfs))  # remove repeated items

records = SeqIO.parse(smorfs, 'fasta')
for record in records:
    smorf_name = record.description
    all_smorfs.append(smorf_name)

    if smorf_name in overlapped_smorfs:
        # print(f"found match {smorf_name}")
        current_overlappedsmorfs.append(smorf_name)

current_overlappedsmorfs_dict = {key: sorted_dict[key] for key in current_overlappedsmorfs if key in sorted_dict}
for key, value in current_overlappedsmorfs_dict.items():
    print(f'{key}: {value}')

# Create a DataFrame from the filtered dictionary
data = {'SMORF': [], 'Overlapping Genes': []}
for key, value in current_overlappedsmorfs_dict.items():
    data['SMORF'].append(key)
    data['Overlapping Genes'].append(", ".join(value))  # join the list of genes into a single string

df_filtered = pd.DataFrame(data)

# Save the DataFrame to an Excel file
df_filtered.to_excel(output_file, index=False)

print(f'Table saved to {output_file}')














