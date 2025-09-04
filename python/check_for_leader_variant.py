from receptor_utils import simple_bio_seq as simple
import glob
import os

bam_fasta_dir = 'f:/ricotta/bam_fastas'
# TRBV12-2*01 in-frame leader variant seen in 10x
leader_seq = 'ATGTAAGATGTGCCTTTTGCTTCCTGCCATGATTCTGAGGCCTCCCCAGCCATGTGGAACT'
gene_seq = 'GATGCTGGCATTATCCAGTCACCCAAGCATGAGGTGACAGAAATGGGACAAACAGTGACTCTGAGATGTGAGCCAATTTTTGGCCACAATTTCCTTTTCTGGTACAGAGATACCTTCGTGCAGGGACTGGAATTGCTGAGTTACTTCCGGAGCCGATCTATTATAGATAATGCAGGTATGCCCACAGAGCGATTCTCAGCTGAGAGGCCTGATGGATCATTCTCTACTCTGAAGATCCAGCCTGCAGAGCAGGGGGACTCGGCCGTGTATGTCT'

# for each .fasta file in bam_fasta_dir, check for the gene seq and, if found, store a sample
found_seqs = {}
seen_seqs = []
for fasta_file in glob.glob(f'{bam_fasta_dir}/*.fasta'):
    sample_name = os.path.basename(fasta_file).replace('.fasta', '')
    if sample_name.startswith('trb'):
        continue
    recs = simple.read_fasta(fasta_file)
    for rec_id, seq in recs.items():
        # find all positions of leader_seq in seq
        l_positions = []
        start = 0
        while True:
            start = seq.find(leader_seq, start)
            if start == -1:
                break
            l_positions.append(start)
            seen_seq = seq[start:start+len(leader_seq)+350]
            if seen_seq not in seen_seqs:
                found_seqs[f'{sample_name}:{rec_id}:{start}'] = seen_seq
                seen_seqs.append(seen_seq)
                print(f'{sample_name}:  {rec_id}  {start}')
            start += 1

simple.write_fasta(f'{bam_fasta_dir}/trbv12-2_leader_variants.fasta', found_seqs)

