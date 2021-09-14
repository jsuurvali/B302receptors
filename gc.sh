#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G

cd /home/suurvalj/scratch/genomes
ls *fna | cut -c-13 | while read asm; do grep -v ">" ${asm}*fna | awk -v x="${asm}" 'BEGIN {FS=""; gc = 0 ; at = 0} {if ($1 ~ /[gcGC]/) gc++ ; else if ($1 ~ /[atAT]/) at++} END {print x" "gc/(gc+at)}' >> gc.tsv; done
