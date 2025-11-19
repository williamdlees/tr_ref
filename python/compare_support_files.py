from receptor_utils import simple_bio_seq as simple

new_file = 'ricotta_aug_28/ricotta_master_allele_support.csv'
old_file = 'ricotta_aug_28/sep_09_version/ricotta_master_allele_support.csv'

new_recs = simple.read_csv(new_file)
old_recs = simple.read_csv(old_file)
new_dict = {rec['vdjbase_name']: rec for rec in new_recs}
old_dict = {rec['vdjbase_name']: rec for rec in old_recs}

all_names = set(new_dict.keys()).union(set(old_dict.keys()))

for name in sorted(all_names):
    if name in new_dict and name in old_dict:
        new_seq = new_dict[name]['seq']
        old_seq = old_dict[name]['seq']
        if new_seq != old_seq:
            print(f'{name} sequence differs')
    elif name in new_dict:
        print(f'{name} only in new file')
    else:
        print(f'{name} only in old file')
