#!/bin/bash
#SBATCH --job-name=proc_PRJNA603748   # Job name
#SBATCH --partition=compute           # Partition or queue name
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=12            # Number of CPU cores per task
#SBATCH --time=6:00:00                # Maximum runtime 
#SBATCH --array=11,13,14,18,20                                                                                                 
#SBATCH --open-mode=truncate                                                                                      
#SBATCH --output=out_files/%a_output.txt
#SBATCH --error=out_files/%a_error.txt 

#Load necessary modules (if needed)
module load anaconda

export PATH=~/cellranger-9.0.1:$PATH

conda activate sra
cd ~/tcr/PRJNA603748

sra_num=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < SRRs.txt)
echo $sra_num

cd ~/tcr/PRJNA603748/${sra_num}

mkdir -p data
fasterq-dump ./${sra_num}.sra --outdir data --split-files --include-technical -p
cd data
mv ${sra_num}_2.fastq ${sra_num}_S1_L001_R1_001.fastq
mv ${sra_num}_3.fastq ${sra_num}_S1_L001_R2_001.fastq
cd ..
cellranger vdj --id PRJNA1083421_${sra_num} --fastqs ./data --reference ~/cellranger_vdj_ref/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --localcores 12 --chain TR
