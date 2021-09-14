#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G

module load hmmer

cd /home/suurvalj/scratch/

echo "shark started at `date`"
hmmsearch -o GCF_902713615.out --domtblout GCF_902713615.domtblout --nobias --nonull2 --domE 0.05 /home/suurvalj/hmm/NLR-TRIM-CASP.hmm GCF_902713615.protein.fa
date
