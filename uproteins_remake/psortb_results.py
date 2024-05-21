"""filter out psortb results"""

import pandas as pd

psortb_dic = {}
psortb_results = "/home/adriana/uproteins/redo/psortb/new_T3_torfs.txt"
df = pd.read_csv(psortb_results, sep="\t", index_col=False)
df_clean = df.fillna(0)
# print(df_clean.head())
# print(df_clean.columns)

for index, row in df_clean.iterrows():
    smorf_names = row['SeqID']
    # print("smorf names:", smorf_names)
    cyt_memb_prot = row['CMSVM+_Localization']
    # print(smorf_names, cyt_memb_prot)
    cell_w_prot = row['CWSVM+_Localization']
    cytop_prot = row['CytoSVM+_Localization']
    extracel_prot = row['ECSVM+_Localization']
    transmem_hel = row['ModHMM+_Details'] # transmembrane helices
    signal_pept = row['Signal+_Details'] # signal peptide
    final_localization = row['Final_Localization']

    # Check if the key has at least one non-"Unknown" value
    if (cyt_memb_prot != "Unknown" or
        cell_w_prot != "Unknown" or
        cytop_prot != "Unknown" or
        extracel_prot != "Unknown" or
        transmem_hel != "No internal helices found" or
        signal_pept != "No signal peptide detected" or
        final_localization != "Unknown"):

        if smorf_names not in psortb_dic:
            psortb_dic[smorf_names] = []

            if cyt_memb_prot != "Unknown":
                psortb_dic[smorf_names].append(cyt_memb_prot)
            if cell_w_prot != "Unknown":
                psortb_dic[smorf_names].append(cell_w_prot)
            if cytop_prot != "Unknown":
                psortb_dic[smorf_names].append(cytop_prot)
            if extracel_prot != "Unknown":
                psortb_dic[smorf_names].append(extracel_prot)
            if transmem_hel != "No internal helices found":
                psortb_dic[smorf_names].append(transmem_hel)
            if signal_pept != "No signal peptide detected":
                psortb_dic[smorf_names].append(signal_pept)
            if final_localization != "Unknown":
                psortb_dic[smorf_names].append(final_localization)

for key, value in psortb_dic.items():
    print(f"{key}: {value}")

for key, values in psortb_dic.items():
    psortb_dic[key] = ', '.join(values)
df_output = pd.DataFrame.from_dict(psortb_dic, orient='index', columns=['Values'])
print(df_output)
df_output.to_csv('/home/adriana/uproteins/redo/psortb/new_T3_torfs_resumed.csv', index_label='SeqID')
