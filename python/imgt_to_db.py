# Create a 'sumary' csv file for TCR sequences. Populate with F and ORF sequences from the IMGT list

from receptor_utils import simple_bio_seq as simple
from Bio import SeqIO
import random
import base64

label_database = list()


def generate_new_label(label_database):
    i = 0
    while i < 10000000:
        s = random.randint(0, 2**20 - 1).to_bytes(5, 'big')
        l = base64.b32encode(s)[4:].decode()
        if l not in label_database:
            label_database.append(l)
            return l

        if i == 100000:
            print('Warning: label namespace is getting full. Perhaps the database is getting too large?')

    print("Error: namespace is reaching exhaustion. Can't allocate a new label in a reasonable amount of time.")
    exit(0)


def read_imgt_fasta(infile: str, species: str, chains=('IGHV', 'IGHD', 'IGHJ', 'CH'), include_orphon: bool =False):
    recs = SeqIO.parse(infile, 'fasta')
    res = []

    for rec in recs:
        desc = rec.description.split('|')
        name = desc[1]
        chain = name[:3]
        type = name[3]
        func = desc[3]
        if ('-REGION' in desc[4]) and chain in chains and desc[2] == species:
            if include_orphon or '/OR' not in name:
                gene, allele = name.split('*')

                res.append({
                    'label': generate_new_label(label_database),
                    'imgt_name': name,
                    'locus': chain,
                    'type': type,
                    'gene': gene,
                    'allele': allele,
                    'imgt_func': func,
                    'imgt_seq': str(rec.seq).upper().replace('.', ''),
                    'seqs': str(rec.seq).upper().replace('.', ''),
                    'longest_seq': str(rec.seq).upper().replace('.', ''),
                    'longest_seq_gapped': str(rec.seq).upper(),
                })
    return res


def main():
    # Read IMGT reference file
    imgt_file = '../../imgt/imgt.fasta'
    summary_file = 'tcr_db.csv'
    recs = read_imgt_fasta(imgt_file, 'Homo sapiens', chains=['TRA', 'TRB', 'TRD', 'TRG'])
    simple.write_csv(summary_file, recs)


main()
