"""filter out smorfs from old interpro results given that we have new tier classification for smorfs now"""

import pandas as pd


df = pd.read_excel('/home/adriana/uproteins/redo/supp_table_1.xlsx', sheet_name=1)
print(df)
smorfs = '/home/adriana/uproteins/redo/240515_smorfTiers_transcriptome.txt'


smorf_dict = {} # contains list of smorfs and their corresponding tiers

with open(smorfs, 'r') as smorf_file:
    next(smorf_file) # skip header

    for line in smorf_file:
        smorf, tier = line.strip().split('\t')
        print("smorf:", smorf)
        smorf_parts = smorf.split('_')
        print(smorf_parts)

        # for gorfs
        # smorf_renamed = "_".join(smorf_parts[0:3])
        # smorf_renamed = smorf_renamed.replace("__", "_")

        # for torfs
        smorf_renamed = "_".join(smorf_parts[:1] + smorf_parts[2:])
        smorf_renamed = smorf_renamed.split('_')
        smorf_renamed = "_".join(smorf_renamed[0:2])
        print("smorf_renamed:", smorf_renamed)


        if tier != 'T4':
            smorf_dict[smorf_renamed] = tier

print(smorf_dict)

filtered_df = df[df.iloc[:, 0].isin(smorf_dict.keys())]
print(filtered_df)
output_filtered_df = '/home/adriana/uproteins/redo/sup_table_1_trans_redo.xlsx'
filtered_df.to_excel(output_filtered_df, index=False)
print(f"Filtered DataFrame has been saved to {output_filtered_df}")