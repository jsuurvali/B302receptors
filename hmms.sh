#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=2
#SBATCH --mem=120G
#SBATCH --job-name=hmms

module load hmmer

cd /home/suurvalj/scratch/translations
cat listall | while read i; do echo "${i} started at `date`"; hmmsearch -o $i.out --domtblout $i.domtblout --nobias --nonull2 --domE 0.05 ~/hmm/NLR-TRIM-CASP.hmm $i.protein.fa; done
