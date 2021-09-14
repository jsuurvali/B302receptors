#!/bin/bash

rm *n25.tsv *n50.tsv; date; for i in *contigs.distrib.tsv; do Rscript contig-counter.R $i 25 >> contigs.distribs.n25.tsv; Rscript contig-counter.R $i 50 >> contigs.distribs.n50.tsv; done; date
for i in *contigs.strands.distrib.tsv; do Rscript contig-counter.R $i 25 >> contigs.strands.distribs.n25.tsv; Rscript contig-counter.R $i 50 >> contigs.strands.distribs.n50.tsv; done; date


grep -P "[0-9]" *contigs.distrib.tsv | sed 's/\.contigs\.distrib\.tsv:/\t/' | grep -v "ll.NACHT" > contigs.distribs.tsv
grep -P "[0-9]" *contigs.strands.distrib.tsv | sed 's/\.contigs\.strands\.distrib\.tsv:/\t/' | grep -v "ll.NACHT" > contigs.strands.distribs.tsv
rm *distrib.tsv


ls contigs*n25.tsv contigs*n50.tsv | while read i; do sed -e 's/GC/\nGC/g' $i | awk 'NR > 1' | tr ' ' '\t' | sed -r 's/\.contigs(\.strands)?\.distrib\.tsv//' > $i; done
