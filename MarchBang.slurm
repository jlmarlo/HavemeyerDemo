#!/bin/bash -l
#SBATCH -t 12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marlo072@umn.edu
#SBATCH --job-name Annotate
#SBATCH -o ./Jobs/%j.Annotate.out
#SBATCH -e ./Jobs/%j.Annotate.err
#SBATCH -p msismall,msilarge

set -e

conda activate snakeprac
#cd $SLURM_SUBMIT_DIR
snakemake -s BandComp.smk -c 6\
	--configfile config.yaml \
	--profile slurm.go_compare \
	--keep-incomplete \
	--rerun-incomplete \
	--use-singularity
