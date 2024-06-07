"""check between old smorf lists and new smorf lists"""

import pandas as pd

old_smorfs_file = '/home/adriana/uproteins/differential_exp/torfs_scts.txt'
new_smorfs_file = '/home/adriana/uproteins/redo/240515_smorfTiers_transcriptome.txt'

new_smorf_list = []
new_smorf_dict = {} # this will contain tier information
df_new = pd.read_csv(new_smorfs_file, sep='\t')
# print(df)
for index, row in df_new.iterrows():
    new_smorf_name = row['smorf']
    new_tier = row['tier']
    if new_smorf_name not in new_smorf_list:
        new_smorf_list.append(new_smorf_name)
        if new_tier != "T4":
            new_smorf_dict[new_smorf_name] = new_tier
# print(new_smorf_dict)


old_smorfs_list = []
df_old = pd.read_csv(old_smorfs_file, sep='\t')
for index, row in df_old.iterrows():
    old_smorf_name = row['full_entries']
    if old_smorf_name not in old_smorfs_list:
        old_smorfs_list.append(old_smorf_name)
# print(old_smorfs_list)

current_smorfs_list = list(set(old_smorfs_list) & set(new_smorf_list))
# for smorf in current_smorfs_list:
    # print(smorf)

# list that will show old smorfs that dont match onto the new smorfs
difference_smorfs_list = list(set(old_smorfs_list) - set(new_smorf_list))
# for smorf in difference_smorfs_list:
#     print(smorf)

final_dictionary = {} # will gather items from current_smorfs_list and their
#corresponding tiers from new_smorfs_dict values
for smorf in current_smorfs_list:
    if smorf in new_smorf_dict:
        final_dictionary[smorf] = new_smorf_dict[smorf]
print(final_dictionary)


output_file = '/home/adriana/uproteins/redo/current_T1-3_torfs_scts.txt'
with open(output_file, 'w') as file:
    for key, value in final_dictionary.items():
        file.write(f"{key}\t{value}\n")

print(f"Final dictionary has been saved to {output_file}")




