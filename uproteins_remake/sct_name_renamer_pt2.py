"""remove other old unannotated genes that aren't supposed to be there"""

import pandas as pd

smorfs_file = "/home/adriana/uproteins/redo/current_T1-3_torfs_SCTs.txt"
counts_file = "/home/adriana/uproteins/isolates/two_strains/two_strains_counts_filtered_updated.txt"

sct_mapping = pd.read_csv(smorfs_file, sep='\t')
data = pd.read_csv(counts_file, sep='\t')

smorf_df = pd.DataFrame(sct_mapping)
gene_sct_list = smorf_df.iloc[:, 0].tolist()
gene_sct_list_modified = [item.replace("gene-", "") for item in gene_sct_list] #remove gene- portion
smorf_sct_list = smorf_df.iloc[:, 2].tolist()
sct_dict = dict(zip(gene_sct_list_modified, smorf_sct_list)) #dictionary containing the sct as keys and their corresponding smorfs as values

count_df = pd.DataFrame(data)
gene_count_list = count_df.iloc[:, 0].tolist()

unannotated_gene_count_list = []
for i in gene_count_list:
    if not i.startswith('Rv'):
        unannotated_gene_count_list.append(i)

# print(unannotated_gene_count_list)

torf_count_list = []
for i in unannotated_gene_count_list:
    if not i.startswith('uproteins'):
        if not i.startswith('tORFs'):
            if not i.startswith('rna'):
                torf_count_list.append(i)

matches = [] # smorfs from T1-T3 that are in counts data
old_smorfs = [] # list of smorfs to be deleted from counts data because they dont belong to T1-T3 anymore
for smorf in torf_count_list:
    if smorf in smorf_sct_list:
        # print(f"found match {smorf}!")
        matches.append(smorf)
    else:
        # print(f"found old smorf {smorf}")
        old_smorfs.append(smorf)

print(old_smorfs)

print("Original dataframe:", count_df)
rows_to_delete = count_df[count_df['gene'].isin(old_smorfs)].index #delete old smorfs
count_df.drop(rows_to_delete, inplace=True)
print("\nUpdated Dataframe:", count_df)

other_genes_to_delete = []
for i in unannotated_gene_count_list:
    if i.startswith('uproteins'):
        other_genes_to_delete.append(i)
    elif i.startswith('tORFs'):
        other_genes_to_delete.append(i)
    # elif i.startswith('rna'):
    #     other_genes_to_delete.append(i)

rows_to_delete2 = count_df[count_df['gene'].isin(other_genes_to_delete)].index #delete rows that have genes that start with 'uproteins', 'rna, or 'torfs'
count_df.drop(rows_to_delete2, inplace=True)
print("\nUpdated Dataframe 2:", count_df)

outfile = "/home/adriana/uproteins/isolates/two_strains/two_strains_counts_filtered_updated_2.xlsx"
count_df.to_excel(outfile, index=False)
print(f"Updated DataFrame saved to {outfile}")



