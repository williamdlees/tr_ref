# compare reference sets in two directories

from receptor_utils import simple_bio_seq as simple

allele_dir_new = 'ricotta_aug_28'
allele_dir_old = 'ricotta_june_7'
loci = ['TRA', 'TRD', 'TRB', 'TRG']
stypes = ['V', 'D', 'J']

for locus in loci:
    for stype in stypes:
        if (locus == 'TRA' or locus == 'TRG') and stype == 'D':
            continue

        old_set = simple.read_fasta(f'{allele_dir_old}/ricotta_{locus}{stype}.fasta')
        new_set = simple.read_fasta(f'{allele_dir_new}/ricotta_{locus}{stype}.fasta')

        old_names = set(old_set.keys())
        new_names = set(new_set.keys())
        shared_names = old_names.intersection(new_names)
        only_old = old_names - new_names
        only_new = new_names - old_names

        print(f'{locus}{stype}: {len(old_names)} old, {len(new_names)} new, {len(shared_names)} shared')

        for name in sorted(shared_names):
            if old_set[name] != new_set[name]:
                print(f'  {name} differs')
