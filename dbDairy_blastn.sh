#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lang.sun@uconn.edu
#SBATCH -o blast.o
#SBATCH -e blast.e


module load blast/2.7.1 

blastn -query allasv.fasta -db DAIRYdb_v1.1.2_blast.fasta -out all.out -max_hsps 5 -outfmt 6 -max_target_seqs 5 -perc_identity 100 
