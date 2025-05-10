from receptor_utils import simple_bio_seq as simple

loci = ['TRA', 'TRD', 'TRB', 'TRG']


def process_locus(locus):
    recs = simple.read_csv(f'ricotta_master_{locus}_alleles.csv')

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
                        'support': [rec['subject']]
                    }
                else:
                    if rec['subject'] not in alleles[rec['vdjbase_allele']]['support']:
                        alleles[rec['vdjbase_allele']]['support'].append(rec['subject'])

    # Collapse support to a count

    for allele in alleles:
        alleles[allele]['support'] = len(alleles[allele]['support'])

    # write summary to csv

    summary = []

    # write out records sorted by stype and allele name
    for allele in sorted(alleles.keys(), key=lambda x: (alleles[x]['stype'], x)):
        summary.append({
            'vdjbase_allele': allele,
            'ttype': alleles[allele]['stype'],
            'support': alleles[allele]['support'],
            'seq': alleles[allele]['seq'],
            'seq_gapped': alleles[allele]['seq_gapped'],
        })

    simple.write_csv(f'riciotta_master_{locus}_allele_support.csv', summary)

    # write fasta files

    for stype in ['V', 'D', 'J']:
        seqs = {k: alleles[k]['seq'] for k in alleles.keys() if alleles[k]['stype'] == stype}
        if seqs:
            simple.write_fasta(f'ricotta_{locus}{stype}.fasta', seqs)
            novels = len({k for k in seqs.keys() if '_' in k})
            print(f'{locus} {stype}: {len(seqs)} / {novels}')

        if stype == 'V':
            seqs = {k: alleles[k]['seq_gapped'] for k in alleles.keys() if alleles[k]['stype'] == stype}
            simple.write_fasta(f'ricotta_{locus}V_gapped.fasta', seqs)


if __name__ == '__main__':
    print('Non-novel/Novel alleles:')
    for locus in loci:
        process_locus(locus)