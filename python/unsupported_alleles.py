# List alleles with too little support to be included in the reference set
# including IMGT alleles that have no support at all

import base64
from turtle import update
from receptor_utils import simple_bio_seq as simple

allele_dir = 'ricotta_aug_28'
loci = ['TRA', 'TRD', 'TRB', 'TRG']
stypes = ['V', 'D', 'J']

rs_fields = [
    'Total_Positions',
    'Average_Coverage',
    'Mismatched_Positions_Coverage_10_Or_Greater',
    'Matched_Positions_Coverage_10_Or_Greater',
    'Position_Mismatches',
    'Position_Matches',
    'Percent_Accuracy',
    'Positions_With_At_Least_10x_Coverage',
    'Fully_Spanning_Reads',
    'Fully_Spanning_Reads_100%_Match'
]


def process_locus(locus):
    supported_recs = simple.read_csv(f'{allele_dir}/ricotta_master_allele_support.csv')
    supported_alleles = [rec['vdjbase_name'] for rec in supported_recs if rec['locus'] == locus]

    recs = simple.read_csv(f'{allele_dir}/allele_table/master_{locus}_alleles.csv')

    alleles = {}
    for rec in recs:
        if rec['vdjbase_allele'] and len(rec['vdjbase_allele']) > 3 and rec['vdjbase_allele'] not in supported_alleles:
            try:
                fm_reads = float(rec['Fully_Spanning_Reads_100%_Match'])
            except ValueError:
                fm_reads = 0

            if rec['vdjbase_allele'] not in alleles:
                stype = rec['vdjbase_allele'][3]
                alleles[rec['vdjbase_allele']] = {
                    'stype': stype,
                    'seq': rec[f'{stype}-REGION'],
                    'seq_gapped': rec['V-REGION-GAPPED'] if stype == 'V' else '',
                    'support': [rec['subject']],
                    'max_matching_reads': int(fm_reads),
                    'max_matching_subj': rec['subject'],
                }
                for field in rs_fields:
                    alleles[rec['vdjbase_allele']][field] = rec[field]
            else:
                if rec['subject'] not in alleles[rec['vdjbase_allele']]['support']:
                    alleles[rec['vdjbase_allele']]['support'].append(rec['subject'])
                if int(fm_reads) > alleles[rec['vdjbase_allele']]['max_matching_reads']:
                    alleles[rec['vdjbase_allele']]['max_matching_reads'] = int(fm_reads)
                    alleles[rec['vdjbase_allele']]['max_matching_subj'] = rec['subject']
                    for field in rs_fields:
                        alleles[rec['vdjbase_allele']][field] = rec[field]

    # Collapse support to a count

    for allele in alleles:
        alleles[allele]['support'] = len(alleles[allele]['support'])

    # create summary

    summary = []

    for allele in sorted(alleles.keys(), key=lambda x: (alleles[x]['stype'], x)):
        if alleles[allele]['stype'] == 'C':
            continue
        
        if len(allele) < 40:
            vdjbase_name = allele
        else:
            vdjbase_name = allele.split('_')[0] + '_' + base64.b32encode(alleles[allele]['seq'].encode('utf-8')).decode()[-4:]

        sum_rec = {
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
        }

        for field in rs_fields:
            sum_rec[field] = alleles[allele][field]

        summary.append(sum_rec)

    return summary, supported_alleles


def main():
    summary = []
    supported_alleles = []
    for locus in loci:
        loc_summary, loc_supported = process_locus(locus)
        summary.extend(loc_summary)
        supported_alleles.extend(loc_supported)

    # now check IMGT alleles for no support at all

    imgt_seqs = {}

    for locus in loci:
        chains = ['TR' + locus[2] + 'V', 'TR' + locus[2] + 'J']
        if locus in ['TRB', 'TRD']:
            chains.append('TR' + locus[2] + 'D')
        recs = simple.read_imgt_fasta('imgt/imgt_set_full.fasta', ['Homo sapiens'], chains, functional_only=True)
        for chain in chains:
            imgt_seqs.update(recs['Homo sapiens'][chain])

    imgt_alleles = imgt_seqs.keys()
    unsupported_alleles = [a['vdjbase_name'] for a in summary]

    no_support = [a for a in imgt_alleles if a not in supported_alleles and a not in unsupported_alleles]

    for allele in sorted(no_support):
        if allele[3] != 'C':
            summary.append({
                'vdjbase_name': allele,
                'locus': allele[:3],
                'type': allele[3],
                'gene': allele.split('*')[0],
                'allele': allele.split('*')[1],
                'support': 0,
                'seq': imgt_seqs[allele],
                'seq_gapped': '',
                'max_matching_reads': 0,
                'max_matching_subj': '',
            })

    simple.write_csv(f'{allele_dir}/ricotta_master_allele_low_or_no_support.csv', summary)


main()
