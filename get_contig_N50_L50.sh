#!/bin/bash

esearch -db assembly -query '"Actinopterygii"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter] AND ("100000"[ContigN50] : "200000000"[ContigN50]))' | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession,Taxid,assembly-status -block Stat -if Stat@category -equals contig_n50 -or Stat@category -equals contig_L50 -sep ":" -def "NA" -element Stat@category,Stat
