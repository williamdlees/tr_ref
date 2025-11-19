# Add columns to the master allele support file indicating matches to any imgt names of sequences including P, ORF, those that were derived from rearranged VDJs - marked in the IMGT file as (F)

from receptor_utils import simple_bio_seq as simple

master_allele_file = 'ricotta_aug_28/ricotta_master_allele_support.csv'
updated_master_allele_file = 'ricotta_aug_28/ricotta_master_allele_support_with_imgt_notes.csv'
allele_dir = 'ricotta_aug_28/allele_table'

imgt_with_rearranged = 'imgt_ref_2024_01_21'

# read imgt refs into dict indexed by sequence

imgt_refs = {}

for locus in ['TRA', 'TRB', 'TRD', 'TRG']:
    for gtype in ['V', 'D', 'J', 'C']:
        if (locus == 'TRA' or locus == 'TRG') and gtype == 'D':
            continue

        imgt_file = f'{imgt_with_rearranged}/Homo_sapiens_{locus}{gtype}.fasta'
        recs = simple.read_fasta(imgt_file)
        for name, seq in recs.items():
            if seq in imgt_refs:
                imgt_refs[seq] = imgt_refs[seq] + ', ' + name
            else:
                imgt_refs[seq] = name


# Read master allele file and update

master_alleles = simple.read_csv(master_allele_file)
vdjbase_to_imgt = {}

seq_lookups = {}
for locus in ['TRA', 'TRB', 'TRD', 'TRG']:
    for gtype in ['V', 'D', 'J', 'C']:
        if (locus == 'TRA' or locus == 'TRG') and gtype == 'D':
            continue

        seq_lookups[f'{locus}{gtype}'] = list()
        for seq, name in imgt_refs.items():
            if name.startswith(f'{locus}{gtype}'):
                seq_lookups[f'{locus}{gtype}'].append((seq, name))

for rec in master_alleles:
    rec['imgt_all_seqs'] = ''
    rec['imgt_sub_super'] = ''
    rec['imgt_sequence'] = ''

    if rec['vdjbase_name'] in vdjbase_to_imgt:
        rec['imgt_all_seqs'], rec['imgt_sub_super'] = vdjbase_to_imgt[rec['vdjbase_name']]
        continue

    if '_' in rec['vdjbase_name']:
        if rec['seq'] in imgt_refs:
            rec['imgt_all_seqs'] = imgt_refs[rec['seq']]

        else:
            # look for sub/super sequences
            loctype = rec['vdjbase_name'][:4]
            for seq, name in seq_lookups[loctype]:
                if rec['seq'] in seq or seq in rec['seq'] and 'C' not in rec['vdjbase_name'] and len(rec['seq']) < 400:
                    comment = ' (imgt is super-seq)' if rec['seq'] in seq else ' (imgt is sub-seq)'
                    if rec['imgt_sub_super']:
                        rec['imgt_sub_super'] = rec['imgt_sub_super'] + ', ' + name + comment
                    else:
                        rec['imgt_sub_super'] = name + comment
                    rec['imgt_sequence'] = seq

    vdjbase_to_imgt[rec['vdjbase_name']] = (rec['imgt_all_seqs'], rec['imgt_sub_super'])
    continue

# filter out obviopusly too long sequences

master_alleles = [rec for rec in master_alleles if (int(rec['Total_Positions']) < 350 or rec['type'] == 'C')]

# Check each novel to make sure its calls in the master files have a reasonable RSS, ie heptamer and nonamer of the right length
# FIrst categorise every novel allele in the master tables as good or bad, where good means at least one example has a reasonable RSS
good_alleles = []
bad_alleles = []
for locus in ['TRA', 'TRB', 'TRD', 'TRG']:
    master_recs = simple.read_csv(f'{allele_dir}/master_{locus}_alleles.csv')
    for rec in master_recs:
        allele = rec['vdjbase_allele']
        if len(allele) < 6:
            continue
        gtype = allele[3]
        if gtype not in ['V', 'J'] or '_' not in allele:
            continue
        if allele in good_alleles:
            continue
        if len(rec[f'{gtype}-HEPTAMER']) == 7 and len(rec[f'{gtype}-NONAMER']) == 9:
            good_alleles.append(allele)
            if allele in bad_alleles:
                bad_alleles.remove(allele)
        else:
            bad_alleles.append(allele)

# Now filter the master alleles
filtered_master_alleles = []
for rec in master_alleles:
    if rec['vdjbase_name'] in bad_alleles:
        print(f'Removing {rec["vdjbase_name"]}: bad RSS')
        continue
    filtered_master_alleles.append(rec)

master_alleles = filtered_master_alleles

# remove TRAJ17*01_a62-_a63- , which has the 4bp deletion affecting its 3' end. We could fix it up as below, but we wouldn't have read support for the added bases
#for rec in master_alleles:
#    if rec['vdjbase_name'] == 'TRAJ17*01_a62-_a63-':
#        rec['seq'] = rec['seq'] + 'GA'
#        rec['vdjbase_name'] = 'TRAJ17*01_a62g'
#        rec['allele'] = '01_a62g'
#        print(f'Fixed up sequence for {rec["vdjbase_name"]} - see notes in google doc')

master_alleles = [rec for rec in master_alleles if rec['vdjbase_name'] != 'TRAJ17*01_a62-_a63-']
print('Removed TRAJ17*01_a62-_a63- which has a 4bp deletion at the 3 end')

simple.write_csv(updated_master_allele_file, master_alleles)

matches = sorted(list(set([rec['imgt_all_seqs'] for rec in master_alleles if rec['imgt_all_seqs'] and 'C' not in rec['vdjbase_name']])))
print(f'Identified {len(matches)} IMGT matches in the larger IMGT set:')
print('\n'.join(matches))
near_matches = list(set([rec['imgt_sub_super'] for rec in master_alleles if rec['imgt_sub_super'] and 'C' not in rec['vdjbase_name'] and len(rec['seq']) < 400]))
print(f'Identified {len(near_matches)} near IMGT rearranged matches:')
displayed_names = []
for rec in master_alleles:
    if rec['imgt_sub_super'] and rec['imgt_sub_super'] not in displayed_names and 'C' not in rec['vdjbase_name'] and len(rec['seq']) < 400:
        displayed_names.append(rec['imgt_sub_super'])
        print(f"{rec['imgt_sub_super']}\n Our sequence: {rec['seq']}\nIMGT sequence: {rec['imgt_sequence']}")
