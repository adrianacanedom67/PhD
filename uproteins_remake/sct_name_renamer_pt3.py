"""rename other scts"""

import pandas as pd

smorfs_file = "/home/adriana/uproteins/redo/current_T1-3_torfs_SCTs.txt"
counts_file = "/home/adriana/uproteins/isolates/two_strains/two_strains_counts_filtered_updated_2.xlsx"

# Try specifying the encoding
sct_mapping = pd.read_csv(smorfs_file, sep='\t', encoding='latin1')
data = pd.read_excel(counts_file, engine='openpyxl')

smorf_df = pd.DataFrame(sct_mapping)
gene_sct_list = smorf_df.iloc[:, 0].tolist()
gene_sct_list_modified = [item.replace("gene-", "") for item in gene_sct_list] # remove gene- portion
smorf_sct_list = smorf_df.iloc[:, 2].tolist()
sct_dict = dict(zip(gene_sct_list_modified, smorf_sct_list)) # dictionary containing the sct as keys and their corresponding smorfs as values

count_df = pd.DataFrame(data)
gene_count_list = count_df.iloc[:, 0].tolist()

scts_to_rename_list = []
for i in gene_count_list:
    if i in gene_sct_list_modified:
        print(f'found {i} in counts data and is a sct')
        scts_to_rename_list.append(i)
    else:
        print("didn't find any scts")

print(scts_to_rename_list)

matching_keys = {}
for key in scts_to_rename_list:
    if key in sct_dict:
        matching_keys[key] = sct_dict[key]

print(matching_keys)

for idx, cell in enumerate(count_df.iloc[:, 0]):
    if cell in matching_keys:
        count_df.iloc[idx, 0] = matching_keys[cell]
print(count_df)

output_renamed_file = "/home/adriana/uproteins/isolates/two_strains/two_strains_counts_filtered_updated_2_RENAMED.xlsx"
count_df.to_excel(output_renamed_file, index=False, engine='openpyxl')
print(f"successfully saved renamed file to {output_renamed_file}")