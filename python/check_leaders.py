from receptor_utils import simple_bio_seq as simple


imgt_v = {}
for locus in ['TRA', 'TRB', 'TRG', 'TRD']:
    recs = simple.read_csv(f'ricotta/alleles_21_aug_2025/master_{locus}_alleles.csv')
    recs = [r for r in recs if r['gene'][3] == 'V']

    for rec in recs:
        if not rec['L-PART1'] or not rec['L-PART2']:
            print(f'Missing L-PART1 or L-PART2 for {rec["subject"]} {rec["gene"]}')
            continue
        exon1 = rec['L-PART1'] + rec['L-PART2']

        # check for whole number of codons
        if len(exon1) % 3 != 0:
            print(f'Exon 1 for {rec["subject"]}{rec["gene"]} is not a whole number of codons')

        # check for stop codon
        s = simple.translate(exon1)
        if '*' in s:
            print(f'Exon 1 {s} for {rec["subject"]} {rec["gene"]} contains a stop codon')