import os
from receptor_utils import simple_bio_seq as simple

proj = 'PRJNA1083421'


def find_fasta_files():
    res = {}
    current_dir = os.getcwd()

    # Iterate through all subdirectories in the current directory
    for subdir in os.listdir(current_dir):
        subdir_path = os.path.join(current_dir, subdir)

        # Check if it's a directory
        if os.path.isdir(subdir_path) and subdir.startswith('SRR'):
            # Construct the expected file path
            nested_dir = f"{proj}_{subdir}"
            fasta_file_path = os.path.join(subdir_path, nested_dir, "outs", "all_contig.fasta")

            # Check if the file exists
            if not os.path.isfile(fasta_file_path):
                print(f"no contig file for sample {subdir}")
                continue

            seqs = simple.read_fasta(fasta_file_path)
            print(f"{subdir}: {len(seqs)} sequences")

            for k, v in seqs.items():
                res[f"{subdir}_{k}"] = v

    return res


if __name__ == "__main__":
    res = find_fasta_files()
    simple.write_fasta(f"{proj}_all_contigs.fasta", res)
    print(f"{len(res)} sequences written to {proj}_all_contigs.fasta")
