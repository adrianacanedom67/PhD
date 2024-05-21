import sys
import os


from Bio.Blast import NCBIXML

def parse(folder, outdir, score_cutoff=50, evalue_cutoff=0.001):
    """
    This function takes a XML file, parses it, and saves the sequences that are similar to the smORFs into a dictionary.
    :param blast_all:
    :param score_cutoff:
    :param evalue_cutoff:
    :return:
    """

    files = os.listdir(folder)
    hits = {}

    for file in files:
        if file.endswith('xml'):
            handler = open(f'{folder}/{file}')
            blast_records = NCBIXML.parse(handler)


            for record in blast_records:
                for hit in record.alignments:
                    for hsp in hit.hsps:
                        # print(vars(record))
                        # print(vars(hsp))  # these are the names inside of each variable on xml file
                        # print(hit.hit_def)
                        taxa = hit.hit_def.split("/")
                        # print(taxa)
                        species = taxa[-1]
                        # print(species)

                        sbjct_sequences = hsp.sbjct
                        # print(sbjct_sequences)

                        score = float(hsp.bits)
                        evalue = float(hsp.expect)
                        if score >= score_cutoff and evalue <= evalue_cutoff:
                            if record.query not in hits:
                                hits[record.query] = []
                                # print(record.query)
                                hits[record.query].append(f'>{record.query}\n{hsp.query}\n')
                            hits[record.query].append(f'>{species}\n{sbjct_sequences}\n')
                            #this will create a dictionary containing all the new smorfs' (other bacteria)
                            #and to which mtb_smorfs they aligned, and since our input files are xml with the other bacteria's
                            #names, the output files will be all the other bacteria's alignments to each smorf
                            #so there'll be a single output file for each mtb_smorf
                            # print(hits[record.query])

    for smorf in hits:
        # print(smorf)
        with open(f'{outdir}/{smorf}.fasta', 'w') as outfile:
            # print(hits[smorf])
            outfile.writelines(hits[smorf])



if __name__ == '__main__':
    parse(folder='/home/farminfo/ACanedo/blastp', outdir='/home/farminfo/ACanedo/smorf_hits')
    # parse(folder=sys.argv[1], outdir=sys.argv[2])