from receptor_utils import simple_bio_seq as simple

# read all IMGT Vs into a dict

imgt_v = {}
for locus in ['TRAV', 'TRBV', 'TRGV', 'TRDV']:
    imgt_v.update(simple.read_fasta(f"imgt/Homo_sapiens_{locus}.fasta"))

ref_v = simple.read_fasta("ricotta_june_7/ricotta_TRV.fasta")

imgt_genes = sorted(list(set([x.split('*')[0] for x in imgt_v.keys()])))
ref_genes = sorted(list(set([x.split('*')[0] for x in ref_v.keys()])))

missing_genes = [x for x in imgt_genes if x not in ref_genes]

for x in missing_genes:
    print(x)
