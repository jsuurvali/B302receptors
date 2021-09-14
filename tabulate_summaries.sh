#!/bin/bash

# converts a set of files containing summaries from nlme::gls output in R to useful tables 

grep -P "(Model|scale\(Genome_len\).*[0-9] .*[0-9] .*[0-9] .*[0-9])" model_summaries.tsv | sed -re 's/  Model: scale\(//' -e 's/\) \~.*//' -e 's/scale\(//' | awk '{if (NR %2 == 1) print ">"$0 ; else print $0}' | tr -d '()' | awk 'BEGIN {RS=">"; OFS="\t"} {print $1,$2,$3,$4,$5,$6}' | sed -e '1d' -e 's/dfdat/all/' -e 's/dat//' | tr -d ')' > model_coefficients.tsv
grep -P "(Model|scale\(lat_c\).*[0-9] .*[0-9] .*[0-9] .*[0-9])" model_summaries.tsv | sed -re 's/  Model: scale\(//' -e 's/\) \~.*//' -e 's/scale\(//' | awk '{if (NR %2 == 1) print ">"$0 ; else print $0}' | tr -d '()' | awk 'BEGIN {RS=">"; OFS="\t"} {print $1,$2,$3,$4,$5,$6}' | sed -e '1d' -e 's/dfdat/all/' -e 's/dat//' | tr -d ')' >> model_coefficients.tsv
grep -P "(Model|scale\(water\).*[0-9] .*[0-9] .*[0-9] .*[0-9])" model_summaries.tsv | sed -re 's/  Model: scale\(//' -e 's/\) \~.*//' -e 's/scale\(//' | awk '{if (NR %2 == 1) print ">"$0 ; else print $0}' | tr -d '()' | awk 'BEGIN {RS=">"; OFS="\t"} {print $1,$2,$3,$4,$5,$6}' | sed -e '1d' -e 's/dfdat/all/' -e 's/dat//' | tr -d ')' >> model_coefficients.tsv


grep -P "(ntercept|Model)" model_summaries.tsv | sed -re 's/scale\(|  Model\:|[\,\~].*//g' | awk '{if (NR %2 == 1) print ">"$0 ; else print $0}'  | tr -d '()' | awk 'BEGIN {RS=">"; OFS="\t"} NR > 1 {print $1,$3,$4,$5,$6}' > model_intercepts.tsv


grep -A3 Model model_summaries.tsv | grep -Pv "(Data|Lik|\-\-)" | sed -re 's/  Model: scale\(//' -e 's/\) \~.*//' -e 's/scale\(//' -e 's/dfdat/all/' -e 's/dat//' | awk '{if (NR %2 == 1) print ">"$0 ; else print $0}' | awk 'BEGIN {RS=">"; OFS="\t"} NR > 1 {print $1,$2,$3,$4}' > model_likelihoods.tsv
grep -P -B2 "(AIC|tandard error)" model_summaries.tsv  | grep -Pv "(Data|^$|\-\-|AIC)" | sed -re 's/Residual standard error\: //' -e 's/  Model: scale\(/>/' -e 's/\).*//' | awk 'BEGIN {RS=">"; OFS="\t"} NR > 1 {print $1,$2,$3,$4,$5,$6,$7}' > model_residuals.tsv
