#!/bin/bash
#SBATCH --job-name=SF3B1_BTK_bc1
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peng.h@wehi.edu.au
#SBATCH --output=find_barcode1.out

module load anaconda3
conda activate /stornext/Genomics/data/CLL_venetoclax/tools/flames_env
module load gcc/11.2.0

for line in HD11 HD12 HD14 CLL63v CLL63i
do echo ${line} ========================
/stornext/Genomics/data/CLL_venetoclax/tools/FLAMES/src/bin/match_cell_barcode \
/vast/scratch/users/peng.h/rachseq/guppy_sup_output/fastq/${line} \
/vast/scratch/users/peng.h/rachseq/guppy_sup_output/fastq/${line}_stat.csv \
/vast/scratch/users/peng.h/rachseq/demultiplex/${line}.fq.gz \
/stornext/Genomics/data/CLL_venetoclax/RaCHseq_paper/barcode/${line}.tsv 2 > /vast/scratch/users/peng.h/rachseq/demultiplex/${line}_demultiplex.log
done
