#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=trnsl

module load emboss

cd /home/suurvalj/scratch/big_genomes/GCA_002915635.3_AmbMex60DD_genomic_big 
ls *.fa | sed 's/\.fa//' | while read i; do echo $i; transeq -sequence $i.fa -outseq $i.protein.fa -frame 6 -clean; done

cd /home/suurvalj/scratch/big_genomes/GCA_016271365.1_neoFor_v3_genomic_big
ls *.fa | sed 's/\.fa//' | while read i; do echo $i; transeq -sequence $i.fa -outseq $i.protein.fa -frame 6 -clean; done
