from Bio.Blast import NCBIXML

result_handle = open('/home/adriana/tcc/smorfs/sequences/renamed_blastp/blastp_finale.xml')

def parse(blast_all):
    blast_records = NCBIXML.parse(result_handle)
    # records = []
    print(result_handle)

    for record in blast_records:
        print(record.query)

        for hit in record.alignments:
            for hsp in hit.hsps:  # matching region between the query and the database hit
                score = float(hsp.score)
                evalue = float(hsp.expect)
                if score >= blast_all.score and evalue <= blast_all.eValue:
                    if record.query not in blast_all.taxa:
                        blast_all.taxa[record.query] = []



