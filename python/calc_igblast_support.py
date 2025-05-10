from receptor_utils import simple_bio_seq as simple

recs = simple.read_csv('PRJNA1083421_all_contigs.tsv', delimiter='\t')

supported_alleles = {}

for rec in recs:
    for stype in ['v', 'd', 'j']:
        if ',' not in rec[f'{stype}_call'] and '100.' in rec[f'{stype}_identity']:
            allele = rec[f'{stype}_call']
            locus = rec['locus']
            if allele not in supported_alleles:
                supported_alleles[allele] = {
                    'locus': locus,
                    'stype': stype.upper(),
                    'support': [rec['sequence_id'].split('_')[0]]
                }
            else:
                sample = rec['sequence_id'].split('_')[0]
                if sample not in supported_alleles[allele]['support']:
                    supported_alleles[allele]['support'].append(sample)

# Collapse support to a count

for allele in supported_alleles:
    supported_alleles[allele]['support'] = len(supported_alleles[allele]['support'])

# write summary to csv

summary = []

# write out records sorted by stype and allele name
for allele in sorted(supported_alleles.keys(), key=lambda x: (supported_alleles[x]['locus'], supported_alleles[x]['stype'], x)):
    summary.append({
        'locus': supported_alleles[allele]['locus'],
        'vdjbase_allele': allele,
        'ttype': supported_alleles[allele]['stype'],
        'support': supported_alleles[allele]['support'],
    })

simple.write_csv(f'PRJNA1083421_all_contigs.csv', summary)

for locus in ['TRA', 'TRD', 'TRB', 'TRG']:
    for stype in ['V', 'D', 'J']:
        seqs = [k for k in supported_alleles.keys() if supported_alleles[k]['stype'] == stype and supported_alleles[k]['locus'] == locus]
        if seqs:
            novels = len([k for k in seqs if '_' in k])
            print(f'{locus} {stype}: {len(seqs)} / {novels}')
