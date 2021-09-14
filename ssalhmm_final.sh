#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --account=def-coling
#SBATCH --cpus-per-task=2
#SBATCH --mem=80G

module load hmmer

# module load emboss

cd /home/suurvalj/scratch/

echo "salmon started at `date` , without looking for zf-FISNACHT"

# wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/237/065/GCA_905237065.2_Ssal_v3.1/GCA_905237065.2_Ssal_v3.1_genomic.fna.gz
# gunzip GCA_905237065.2_Ssal_v3.1_genomic.fna.gz
# transeq -sequence GCA_905237065.2_Ssal_v3.1_genomic.fna -outseq GCA_905237065.protein.fa -clean -frame 6
hmmsearch -o GCA_905237065.out --domtblout GCA_905237065.domtblout --nobias --nonull2 --domE 0.05 /home/suurvalj/hmm/NLR-TRIM-CASP_no-zfFISNACHT.hmm GCA_905237065.protein.fa

date
