# Create reference set from genomic allele table
# Also create a 'support' file listing all alleles with their support counts
# Update the db with the alleles

from receptor_utils import simple_bio_seq as simple
import random
import base64

allele_dir = 'ricotta_june_2'
tcr_db_name_in = 'db/2025_06_04_imgt_start/tcr_db.csv'
tcr_db_name_out = 'db/2025_06_04/tcr_db.csv'

label_database = []


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


loci = ['TRA', 'TRD', 'TRB', 'TRG']

substitution_dict = {
    "TRAV1-1*01_c170g": "TRAV1-1*02",
    "TRAV12-2*01_g86t_t164c": "TRAV12-2*03",
    "TRAV12-2*01_t164c": "TRAV12-2*02",
    "TRAV14/DV4*01_c306t": "TRAV14/DV4*03",
    "TRAV21*01_g15a": "TRAV21*02",
    "TRAV35*01_a135c": "TRAV35*02",
    "TRAV36/DV7*05_c54t_a301g": "TRAV36/DV7*04",
    "TRAV38-1*01_a89g_g129a": "TRAV38-1*03",
    "TRAV6*01_c87t": "TRAV6*02",
    "TRAV6*01_c87t_a271g": "TRAV6*04",
    "TRAV6*01_c87t_c165t": "TRAV6*03",
    "TRAV8-4*01_a45g_t172a_c176g": "TRAV8-4*03",
    "TRAV8-7*01_324caggag329": "TRAV8-7*01",
    "TRAV8-7*02_324caggag329": "TRAV8-7*02",
    "TRAV9-2*02_a246g": "TRAV9-2*03",
    "TRDV2*01_g49a": "TRDV2*02",
    "TRBV10-3*01_t114c_g156a": "TRBV10-3*03",
    "TRBV11-3*01_t128g": "TRBV11-3*04",
    "TRBV14*01_g204a": "TRBV14*02",
    "TRBV19*01_g171c": "TRBV19*03",
    "TRBV2*01_g65a": "TRBV2*02",
    "TRBV20-1*01_t28a": "TRBV20-1*02",
    "TRBV20-1*01_t28a_c142a": "TRBV20-1*05",
    "TRGV4*02_27d6": "TRGV2*02",
}

# also don't coalesce Ds and Js.

def process_locus(locus):
    recs = simple.read_csv(f'{allele_dir}/allele_table/master_{locus}_alleles.csv')

    alleles = {}
    for rec in recs:
        if rec['vdjbase_allele'] and len(rec['vdjbase_allele']) > 3:
            try:
                fs_reads = float(rec['Fully_Spanning_Reads'])
                fm_reads = float(rec['Fully_Spanning_Reads_100%_Match'])
            except ValueError:
                continue

            if fs_reads > 7 and fm_reads/fs_reads >= 0.8:
                if rec['vdjbase_allele'] not in alleles:
                    stype = rec['vdjbase_allele'][3]
                    alleles[rec['vdjbase_allele']] = {
                        'stype': stype,
                        'seq': rec[f'{stype}-REGION'],
                        'seq_gapped': rec['V-REGION-GAPPED'] if stype == 'V' else '',
                        'support': [rec['subject']],
                        'max_matching_reads': int(float(rec['Fully_Spanning_Reads_100%_Match'])),
                        'max_matching_subj': rec['subject'],
                    }
                else:
                    if rec['subject'] not in alleles[rec['vdjbase_allele']]['support']:
                        alleles[rec['vdjbase_allele']]['support'].append(rec['subject'])
                    if int(float(rec['Fully_Spanning_Reads_100%_Match'])) > alleles[rec['vdjbase_allele']]['max_matching_reads']:
                        alleles[rec['vdjbase_allele']]['max_matching_reads'] = int(float(rec['Fully_Spanning_Reads_100%_Match']))
                        alleles[rec['vdjbase_allele']]['max_matching_subj'] = rec['subject']

    # Collapse support to a count

    for allele in alleles:
        alleles[allele]['support'] = len(alleles[allele]['support'])

    # write summary to csv

    summary = []

    for allele in sorted(alleles.keys(), key=lambda x: (alleles[x]['stype'], x)):
        if len(allele) < 40:
            vdjbase_name = allele
        else:
            vdjbase_name = allele.split('_')[0] + '_' + base64.b32encode(alleles[allele]['seq'].encode('utf-8')).decode()[-4:]

        summary.append({
            'vdjbase_name': vdjbase_name,
            'locus': locus,
            'type': alleles[allele]['stype'],
            'gene': allele.split('*')[0],
            'allele': allele.split('*')[1],
            'support': alleles[allele]['support'],
            'seq': alleles[allele]['seq'],
            'seq_gapped': alleles[allele]['seq_gapped'],
            'max_matching_reads': alleles[allele]['max_matching_reads'],
            'max_matching_subj': alleles[allele]['max_matching_subj'],
        })

    for stype in ['V', 'D', 'J']:
        seqs = {k: alleles[k]['seq'] for k in alleles.keys() if alleles[k]['stype'] == stype}
        if seqs:
            simple.write_fasta(f'{allele_dir}/ricotta_{locus}{stype}.fasta', seqs)

        if stype == 'V':
            seqs = {k: alleles[k]['seq_gapped'] for k in alleles.keys() if alleles[k]['stype'] == stype}
            simple.write_fasta(f'{allele_dir}/ricotta_{locus}V_gapped.fasta', seqs)

    return summary


def main():
    summary = []
    for locus in loci:
        summary.extend(process_locus(locus))

    simple.write_csv(f'{allele_dir}/ricotta_master_allele_support.csv', summary)

    # pan-locus sets

    for stype in ['V', 'D', 'J']:
        seqs = {rec['vdjbase_name']: rec['seq'] for rec in summary if rec['type'] == stype}
        if seqs:
            simple.write_fasta(f'{allele_dir}/ricotta_TR{stype}.fasta', seqs)

        if stype == 'V':
            seqs = {rec['vdjbase_name']: rec['seq_gapped'] for rec in summary if rec['type'] == 'V'}
            simple.write_fasta(f'{allele_dir}/ricotta_TRV_gapped.fasta', seqs)

    # Update the database with the new alleles
    tcr_db = simple.read_csv(tcr_db_name_in)
    label_database = [rec['label'] for rec in tcr_db]
    tcr_db = {rec['longest_seq']: rec for rec in tcr_db}

    for rec in summary:
        found = False
        if rec['seq'] in tcr_db:
            tcr_db[rec['seq']]['vdjbase_name'] = rec['vdjbase_name']
            tcr_db[rec['seq']]['supporting_subjs'] = rec['support']
            tcr_db[rec['seq']]['max_matching_reads'] = rec['max_matching_reads']
            tcr_db[rec['seq']]['max_matching_subj'] = rec['max_matching_subj']
            found = True

        # is this a super-sequence?
        if not found:
            for seq in list(tcr_db.keys()):
                if seq in rec['seq']:
                    tcr_db[seq]['vdjbase_name'] = rec['vdjbase_name']
                    tcr_db[seq]['max_matching_reads'] = rec['max_matching_reads']
                    tcr_db[seq]['max_matching_subj'] = rec['max_matching_subj']
                    tcr_db[seq]['supporting_subjs'] = rec['support']
                    tcr_db[seq]['seqs'] = ','.join(tcr_db[seq]['seqs'].split(',') + [rec['seq']])
                    tcr_db[seq]['longest_seq'] = rec['seq']
                    if 'notes' not in tcr_db[seq]:
                        tcr_db[seq]['notes'] = ''
                    tcr_db[seq]['notes'] += f'genomic allele is a super-sequence'
                    tcr_db[rec['seq']] = tcr_db[seq]
                    del tcr_db[seq]
                    found = True
                    break

        # is this a sub-sequence?
        if not found:
            for seq in tcr_db:
                if rec['seq'] in seq:
                    tcr_db[seq]['vdjbase_name'] = rec['vdjbase_name']
                    tcr_db[seq]['max_matching_reads'] = rec['max_matching_reads']
                    tcr_db[seq]['max_matching_subj'] = rec['max_matching_subj']
                    tcr_db[seq]['supporting_subjs'] = rec['support']
                    tcr_db[seq]['seqs'] = ','.join(tcr_db[seq]['seqs'].split(',') + [rec['seq']])
                    if 'notes' not in tcr_db[seq]:
                        tcr_db[seq]['notes'] = ''
                    tcr_db[seq]['notes'] += f'genomic allele is a sub-sequence'
                    found = True
                    break

        if not found:
            tcr_db[rec['seq']] = {
                'label': generate_new_label(label_database),
                'vdjbase_name': rec['vdjbase_name'],
                'locus': rec['locus'],
                'type': rec['type'],
                'gene': rec['gene'],
                'allele': rec['allele'],
                'imgt_func': '',
                'imgt_seq': '',
                'seqs': rec['seq'],
                'longest_seq': rec['seq'],
                'longest_seq_gapped': rec['seq_gapped'],
                'supporting_subjs': rec['support'],
                'max_matching_reads': rec['max_matching_reads'],
                'max_matching_subj': rec['max_matching_subj'],
            }

    # Write the updated database
    simple.write_csv(tcr_db_name_out, list(tcr_db.values()), scan_all=True)


main()
