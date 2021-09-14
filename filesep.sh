#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=fsep

cd /home/suurvalj/scratch/big_genomes/neofor
grep ">" GCA_016271365.1_neoFor_v3_genomic_big.fna | tr -d ">" | while read i; do awk -v x="${i}" 'BEGIN {RS=">" ; FS="\n"} $1 == x {print ">"$0}' GCA_016271365.1_neoFor_v3_genomic_big.fna > $i.fa; done
