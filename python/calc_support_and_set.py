from receptor_utils import simple_bio_seq as simple

loci = ['TRA', 'TRD', 'TRB', 'TRG']


def process_locus(locus):
    print(f'\n{locus}...')
    recs = simple.read_csv(f'master_{locus}_alleles.csv')

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

    simple.write_csv(f'ricotta_master_{locus}_allele_support.csv', summary)

    # write fasta files

    print("Number of alleles in genomic data non-novel/novel:")

    for stype in ['V', 'D', 'J']:
        seqs = {k: alleles[k]['seq'] for k in alleles.keys() if alleles[k]['stype'] == stype}
        if seqs:
            simple.write_fasta(f'ricotta_{locus}{stype}.fasta', seqs)
            novels = len({k for k in seqs.keys() if '_' in k})
            print(f'{locus} {stype}: {len(seqs)-novels} / {novels}')

        if stype == 'V':
            seqs = {k: alleles[k]['seq_gapped'] for k in alleles.keys() if alleles[k]['stype'] == stype}
            simple.write_fasta(f'ricotta_{locus}V_gapped.fasta', seqs)

    # imgt comparison
    print("IMGT allele names found/not found in genomic data:")
    imgt_recs = simple.read_fasta('../imgt/imgt_TR.fasta')
    imgt_recs = {k: v for k, v in imgt_recs.items() if k.startswith(locus)}
    allele_names = list(alleles.keys())

    for type in ['V', 'D', 'J']:
        if type == 'D' and locus not in ('TRB', 'TRD'):
            continue
        imgt_found = 0
        imgt_notfound = 0

        for imgt_name in imgt_recs.keys():
            if imgt_name.startswith(f'{locus}{type}'):
                if imgt_name in allele_names:
                    imgt_found += 1
                else:
                    imgt_notfound += 1
        print(f'{locus} {type}: {imgt_found}/{imgt_notfound}')

    # sequence comparison
    print("IMGT sequencea matched/partially matched/not matched in genomically derived sequences:")
    res = []
    for type in ['V', 'D', 'J']:
        if type == 'D' and locus not in ('TRB', 'TRD'):
            continue
        imgt_seqs = {k: v for k, v in imgt_recs.items() if k.startswith(f'{locus}{type}')}
        len_all = len(imgt_seqs)

        if type == 'V':
            imgt_seqs = {k: v for k, v in imgt_seqs.items() if '*' not in simple.translate(v)}
            print(f'Filtering out IMGT V alleles with stop codons: {len_all} -> {len(imgt_seqs)}')

        allele_seqs = {k: alleles[k]['seq'] for k in alleles.keys() if alleles[k]['stype'] == type}
        exact = 0
        partial = 0
        nonmatch = 0
        
        for imgt_name, imgt_seq in imgt_seqs.items():
            matched_name = ''
            rec = {'imgt_name': imgt_name, 'imgt_seq': imgt_seq, 'match': '', 'genomic_name': '', 'genomic_seq': ''}
            if imgt_seq in allele_seqs.values():
                matched_name = [k for k, v in allele_seqs.items() if v == imgt_seq][0]
                rec['match'] = 'exact'
                exact += 1
            else:
                for allele_name, allele_seq in allele_seqs.items():
                    if imgt_seq in allele_seq or allele_seq in imgt_seq:
                        matched_name = allele_name
                        partial += 1
                        rec['match'] = 'partial'
                        break

            if not matched_name:
                rec['match'] = 'nonmatch'
                nonmatch += 1

            rec['genomic_name'] = matched_name
            if matched_name:
                rec['genomic_seq'] = allele_seqs[matched_name]

            res.append(rec)

        print(f'{locus} {type}: {exact}/{partial}/{nonmatch}')
    simple.write_csv(f'ricotta_{locus}_allele_imgt_seq_comparison.csv', res)


if __name__ == '__main__':
    print('Total/Novel alleles:')
    for locus in loci:
        process_locus(locus)
