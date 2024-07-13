"""check which gorfs/torfs are equivalent based on their sequences, one fasta at a time"""

from Bio import SeqIO


smorf_hits = "/home/adriana/uproteins/redo/phylogenetic_analysis/old_xml/tblastn/gorfs_T5_2.txt"
fasta_file = "/home/adriana/tcc/smorfs/sequences/gorfs_T5.fasta"
torfs_file = "/home/adriana/tcc/smorfs/sequences/torfs_T1-5.fasta"


def hit_gatherer(smorf_hits):
    hits_list = []
    with open(smorf_hits, 'r') as txt_file:
        for line in txt_file:
            hits_list.append(line.strip())
    return hits_list


def check_smorf_in_hits(fasta_file, hits_list):
    records = SeqIO.parse(fasta_file, 'fasta')
    hits_dict = {}
    for record in records:
        smorf = record.description
        sequence = str(record.seq)
        found = False
        for hit in hits_list:
            if hit in smorf:
                found = True
                break
        if found:
            # print(f'{smorf} found')
            hits_dict[smorf] = sequence
    # print(hits_dict)
    return hits_dict

def sequence_checker(torfs_file, hits_dict):
    records = SeqIO.parse(torfs_file, 'fasta')
    torfs_dict = {}
    for record in records:
        torf = record.description
        sequence = str(record.seq)
        torfs_dict[sequence] = torf
    same_sequences = []
    for smorf, hit_sequence in hits_dict.items():
        if hit_sequence in torfs_dict:
            same_sequences.append(smorf)

    print(same_sequences)
    return same_sequences



hits_list = hit_gatherer(smorf_hits)
hits_dict = check_smorf_in_hits(fasta_file, hits_list)
same_sequences = sequence_checker(torfs_file, hits_dict)