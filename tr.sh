#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

module load emboss

# This script reads in all *.fa files in a folder and translates them in six reading frames
# All stop codons are replaced with the character X
# The output files are named *.protein.fa
# The "#SBATCH" and "module load" lines are for running the script on SLURM

cd /home/suurvalj/scratch/big_genomes      # replace this with actual location of the files
ls *.fa | sed 's/\.fa//' | while read i; do echo $i; transeq -sequence $i.fa -outseq $i.protein.fa -frame 6 -clean; done
