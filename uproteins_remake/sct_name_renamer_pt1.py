"""rename SCTs according to smorf names"""

import pandas as pd
import sys

smorfs_file = r'C:\Users\LENOVO\Desktop\PhD\uproteins_paper\clinical_isolates\current_T1-3_torfs_SCTs.txt'
counts_file = r'C:\Users\LENOVO\Desktop\PhD\uproteins_paper\clinical_isolates\two_strains_counts_filtered.txt'

data = pd.read_csv(counts_file, sep='\t') # read rna-seq counts file
sct_mapping = pd.read_csv(smorfs_file, sep='\t') # read smorfs file

count_df = pd.DataFrame(data)
gene_list = count_df.iloc[:, 0].tolist()

smorf_df = pd.DataFrame(sct_mapping)
sct_list = smorf_df.iloc[:, 0].tolist()
sct_modified_list = [item.replace("gene-", "") for item in sct_list] # remove "gene-" portion of list
smorf_list = smorf_df.iloc[:, 2].tolist()

# dictionary containing all the Rvs as keys and the corresponding smorfs as values
sct_dict = dict(zip(sct_modified_list, smorf_list))
# print(sct_dict)

to_delete = [] # list containing those genes that do not belong to T1-T3 anymore and must be deleted from the counts data

for i in gene_list:
    if not i.startswith('Rv'): # gather which unnanotated genes are already in this table
        if not i in smorf_list:
            to_delete.append(i)

# delete genes from counts table that don't belong to T1-T3 anymore
count_df_filtered = count_df[~count_df.iloc[:, 0].isin(to_delete)]
# print(count_df_filtered)
# print(count_df)
counts_df_filtered_output = r'C:\Users\LENOVO\Desktop\PhD\uproteins_paper\clinical_isolates\two_strains_counts_filtered_updated.txt'
print(f"successfully saved file to {counts_df_filtered_output}")






