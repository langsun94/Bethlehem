#!/bin/bash
#SBATCH --job-name=sourcetracker
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lang.sun@uconn.edu
#SBATCH -o %x_%j.eut
#SBATCH -e %x_%j.err

module load R/3.5.1
module load qiime/1.9.1

filter_otus_from_otu_table.py -i ps.fungi.0.biom -o filtered_otu_table.biom -s 1

biom convert -i filtered_otu_table.biom -o table.from_biom.txt --to-tsv

R --slave --vanilla --args -i table.from_biom.txt -m rawmilk.txt -o rawmilk < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

R --slave --vanilla --args -i table.from_biom.txt -m unripened.txt -o unripened < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

R --slave --vanilla --args -i table.from_biom.txt -m ripened.txt -o ripened < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

R --slave --vanilla --args -i table.from_biom.txt -m mixed.txt -o mixed < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

R --slave --vanilla --args -i table.from_biom.txt -m curd.txt -o curd < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

R --slave --vanilla --args -i table.from_biom.txt -m cheese.txt -o cheese < $SOURCETRACKER_PATH/sourcetracker_for_qiime.r

