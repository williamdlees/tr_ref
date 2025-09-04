# Create reference set from genomic allele table
# Also create a 'support' file listing all alleles with their support counts
# Update the db with the alleles

from receptor_utils import simple_bio_seq as simple
import random
import base64

allele_dir = 'ricotta_aug_28'
tcr_db_name_in = 'db/2025_06_07/tcr_db.csv'
tcr_db_name_out = 'db/2025_08_28/tcr_db.csv'

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

        # is this a V super-sequence?
        if not found and rec['type'] == 'V':
            for seq in list(tcr_db.keys()):
                if seq in rec['seq'] and tcr_db[seq]['type'] == 'V':
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
                    if seq == 'ACCCAGTCGGTGACCCAGCTTGATGGCCACATCACTGTCTCTGAAGAAGCCCCTCTGGAACTGAAGTGCAACTATTCCTATAGTGGAGTTCCTTCTCTCTTCTGGTATGTCCAATACTCTAGCCAAAGCCTCCAGCTTCTCCTCAAAGACCTAACAGAGGCCACCCAGGTTAAAGGCATCAGAGGTTTTGAGGCTGAATTTAAGAAGAGCGAAACCTCCTTCTACCTGAGGAAACCATCAACCCATGTGAGTGATGCTGCTGAGTACTTCTGTGCTGTGGGTGACAGGAG':
                        breakpoint()
                    del tcr_db[seq]
                    found = True
                    break

        # is this a V sub-sequence?
        # if so make a new record
        if not found and rec['type'] == 'V':
            for seq in tcr_db:
                if rec['seq'] in seq:
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
                        'longest_seq_gapped': rec['seq_gapped'] if rec['seq_gapped'] else rec['seq'],
                        'supporting_subjs': rec['support'],
                        'max_matching_reads': rec['max_matching_reads'],
                        'max_matching_subj': rec['max_matching_subj'],
                        'notes': 'sub-sequence of ' + tcr_db[seq]['label'],
                    }
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
                'longest_seq_gapped': rec['seq_gapped'] if rec['seq_gapped'] else rec['seq'],
                'supporting_subjs': rec['support'],
                'max_matching_reads': rec['max_matching_reads'],
                'max_matching_subj': rec['max_matching_subj'],
            }

    # Write the updated database
    simple.write_csv(tcr_db_name_out, list(tcr_db.values()), scan_all=True)

    # pan-locus reference sets

    for stype in ['V', 'D', 'J']:
        seqs = {}
        for rec in tcr_db.values():
            if rec['type'] == stype and 'vdjbase_name' in rec and rec['vdjbase_name']:
                name = rec['imgt_name'] if ('imgt_name' in rec and rec['imgt_name']) else rec['vdjbase_name']
                seqs[name] = rec['longest_seq']
        if seqs:
            seqs = {k: seqs[k] for k in sorted(seqs.keys())}
            simple.write_fasta(f'{allele_dir}/ricotta_TR{stype}.fasta', seqs)

        if stype == 'V':
            seqs = {}
            for rec in tcr_db.values():
                if rec['type'] == 'V' and 'vdjbase_name' in rec and 'vdjbase_name' in rec and rec['vdjbase_name']:
                    name = rec['imgt_name'] if ('imgt_name' in rec and rec['imgt_name']) else rec['vdjbase_name']
                    seqs[name] = rec['longest_seq_gapped']
            seqs = {k: seqs[k] for k in sorted(seqs.keys())}
            simple.write_fasta(f'{allele_dir}/ricotta_TRV_gapped.fasta', seqs)

    # per-locus reference sets
    for locus in loci:  
        for stype in ['V', 'D', 'J']:
            seqs = {}
            for rec in tcr_db.values():
                if rec['locus'] == locus and rec['type'] == stype and 'vdjbase_name' in rec and rec['vdjbase_name']:
                    name = rec['imgt_name'] if ('imgt_name' in rec and rec['imgt_name']) else rec['vdjbase_name']
                    seqs[name] = rec['longest_seq']
            if seqs:
                # sort by name
                seqs = {k: seqs[k] for k in sorted(seqs.keys())}
                simple.write_fasta(f'{allele_dir}/ricotta_{locus}{stype}.fasta', seqs)

            if stype == 'V':
                seqs = {}
                for rec in tcr_db.values():
                    if rec['locus'] == locus and rec['type'] == 'V' and 'vdjbase_name' in rec and rec['vdjbase_name']:
                        name = rec['imgt_name'] if ('imgt_name' in rec and rec['imgt_name']) else rec['vdjbase_name']
                        seqs[name] = rec['longest_seq_gapped']
                seqs = {k: seqs[k] for k in sorted(seqs.keys())}        
                simple.write_fasta(f'{allele_dir}/ricotta_{locus}V_gapped.fasta', seqs)

main()
