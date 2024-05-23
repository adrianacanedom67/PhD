import sys

import pandas as pd
from Bio import SeqIO

def fasta_maker(file, sequences):
    df = pd.read_excel(file)

    protein_names_list = df.iloc[:, 0].tolist()
    print(protein_names_list)

    protein_sequences = {} # dictionary to store the proteins found with their corresponding sequences

    records = SeqIO.parse(sequences, 'fasta')
    for record in records:
        identifiers = str(record.description)
        sequences = str(record.seq)
        # print("fasta proteins:", identifiers)
        # print("identifiers:", identifiers)
        protein_name_parts = identifiers.split(' ')
        protein_name_in_fasta = protein_name_parts[0]
        # print("protein name split:", protein_name_in_fasta)
        if protein_name_in_fasta in protein_names_list:
            print(f"found protein {protein_name_in_fasta} in list with sequence")
            protein_sequences[protein_name_in_fasta] = sequences

    print(protein_sequences)

    # filename = file.replace(".xlsx", "")
    fasta_output_path = f'/home/adriana/covid/folder/Projeto2/Mock_vs_Infected7rep_notL_Stringency.fasta'
    with open(fasta_output_path, 'w') as fasta_file:
        for protein_name_in_fasta, sequences in protein_sequences.items():
            fasta_file.write(f'>{protein_name_in_fasta}\n')
            fasta_file.write(f'{sequences}\n')



if __name__ == '__main__':
    fasta_maker(file=sys.argv[1],
                sequences='/home/adriana/covid/homosapiens+sarscov2_proteome_UP000005640_2024_04_21.fasta')
