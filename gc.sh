#!/bin/bash

# calculates the GC content of a set of files that end with *fna (genome sequences from the NCBI databases).
# The first 13 characters of each file name are used to label the output.

ls *fna | cut -c-13 | while read asm; do \
  grep -v ">" ${asm}*fna | \
  awk -v x="${asm}" 'BEGIN {FS=""; gc = 0 ; at = 0} {if ($1 ~ /[gcGC]/) gc++ ; else if ($1 ~ /[atAT]/) at++} END {print x" "gc/(gc+at)}' \
  >> gc.tsv
done
