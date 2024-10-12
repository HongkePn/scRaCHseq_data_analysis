#!/bin/bash
#SBATCH --job-name=guppy_sup_basecall
#SBATCH --time=48:00:00
#SBATCH --output=out.out
#SBATCH --error=err.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peng.h@wehi.edu.au
#SBATCH --cpus-per-task=16
#SBATCH --partition=gpuq
#SBATCH --gres=gpu:A30:4
#SBATCH --mem=40G


module purge
module load guppy-gpu/6.2.1
guppy_basecaller \
--input_path /vast/scratch/users/peng.h/rachseq/fast5/fast5_pass \
--save_path /vast/scratch/users/peng.h/rachseq/guppy_sup_output \
-c dna_r9.4.1_450bps_sup_prom.cfg \
--barcode_kits SQK-PCB111-24 \
--recursive \
-x 'cuda:all' \
--num_callers 4 \
--min_qscore 9 \
--compress_fastq \
--chunks_per_runner 128
