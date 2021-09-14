## PRIOR TO RUNNING, PREPARE THE DATA:
## cat domains_red.tsv | while read i; do awk 'BEGIN {OFS="\t"} $1 !~/#/ && $17-$16+1 >= 30 && $13 < 1e-05 && $4 != "zf_FISNACHT" {print $1,$4,$18,$19,$6,$8}' GCF_000002035.domtblout | sed -re 's/^([^\t]+_)([1-3])\t/\1fw\t\2\t/' -e 's/^([^\t]+_)([4-6])\t/\1rv\t\2\t/' -e 's/(zf_B30\.2|PRY|SPRY)/B30\.2/' -e 's/zf-B_box/B-BOX/' -e 's/(zf-C3HC4_?[0-9]?[0-9]?|zf-RING_UBOX)/C3HC4/' -e 's/CARD_2/CARD/' | awk -v x="${i}" '$3 == x && ! ($5 - $4 < 100  && $7 =="inf" && 0.9*$6 > $5-$4) {if ($2 <= 3 ) print $1"\t"($4-1)*3+$2"\t"$5*3-1+$2; else if ($2 >= 4) print $1"\t"$4*3+$2-6"\t"$5*3+$2-4}' | bedtools sort |bedtools merge -d 30 | awk -v x="${i}" 'BEGIN {OFS="\t"} {print $1,$2+1,$3,x}' >> /mnt/e/Desktop/tmp.tsv; done

# next two rows are just there in case one needs to find the length of a chromosome or contig quickly
# library(rentrez)
# entrez_summary(db="nuccore", id = "contig_id")$slen

setwd("E:/Dropbox/TRIMpaper_2021/preprocessed")
library(data.table)

assemblies = read.table("assemblies.txt")[,1] # this is actually just a list of identifiers so that the script would find the files
for (assembly in assemblies) {
  subsets <- c("all.NACHT",
               "all.B302",
               "all.BBOX",
               "all.TRIM.InclSoloBBOX",
               "all.TRIM",
               "FISNA.NACHT",
               "TRIM.fn3",
               "TRIM.B302",
               "NLR.B302",
               "NLR.FIIND.CARD",
               "CARD.NLR",
               "PYD.NLR",
               "PYD.NLR.B302",
               "CARD.NLR.B302")
  
  
# read in the data and sort it by genomic coordinates
  
  asmfilename <- paste(assembly, ".tsv", sep = "")
  dt <- fread(asmfilename)
  
  colnames(dt) <- c("contig", "current.domain.start", "current.domain.end", "current.domain.type")
  dt <- dt[order(contig, current.domain.start)]
  
  # get the contents of the previous and next row
  
  dt[ , previous.contig := shift(contig, type = 'lag')]
  dt[ , previous.domain.start := shift(current.domain.start, type = 'lag')]
  dt[ , previous.domain.end := shift(current.domain.end, type = 'lag')]
  dt[ , previous.domain.type := shift(current.domain.type, type = 'lag')]
  
  dt[ , next.contig := shift(contig, type = 'lead')]
  dt[ , next.domain.start := shift(current.domain.start, type = 'lead')]
  dt[ , next.domain.end := shift(current.domain.end, type = 'lead')]
  dt[ , next.domain.type := shift(current.domain.type, type = 'lead')]
  
  
  # if the previous or next field is not within 100kb, replace contents with NA
  
  dt$previous.contig <- ifelse(dt$contig != dt$previous.contig | dt$current.domain.start - dt$previous.domain.end > 99999, NA, dt$previous.contig)
  dt$previous.domain.start <- ifelse(is.na(dt$previous.contig), NA, dt$previous.domain.start)
  dt$previous.domain.end <- ifelse(is.na(dt$previous.contig), NA, dt$previous.domain.end)
  dt$previous.domain.type <- ifelse(is.na(dt$previous.contig), NA, dt$previous.domain.type)
  
  dt$next.contig <- ifelse(dt$contig != dt$next.contig | dt$next.domain.start - dt$current.domain.end > 99999, NA, dt$next.contig)
  dt$next.domain.start <- ifelse(is.na(dt$next.contig), NA, dt$next.domain.start)
  dt$next.domain.end <- ifelse(is.na(dt$next.contig), NA, dt$next.domain.end)
  dt$next.domain.type <- ifelse(is.na(dt$next.contig), NA, dt$next.domain.type)
  
  
  setcolorder(dt, c(1,8,6,7,4,2,3,12,10,11))
  dt <- dt[, ! c("previous.contig", "next.contig"), with = FALSE]
  
  
  out.alldat <- list()
  
  out.alldat[[1]] <- dt[current.domain.type == "NACHT"]
  out.alldat[[2]] <- dt[current.domain.type == "B30.2"]
  out.alldat[[3]] <- dt[current.domain.type == "B-BOX"]
  out.alldat[[4]] <- dt[current.domain.type == "B-BOX" | (previous.domain.type == "C3HC4" & current.domain.type %in% c("B30.2", "fn3"))]
  out.alldat[[5]] <- dt[(current.domain.type == "B-BOX" & next.domain.type %in% c("B30.2", "fn3")) | (previous.domain.type == "C3HC4" & current.domain.type %in% c("B-BOX", "B30.2", "fn3"))]
  out.alldat[[6]] <- dt[current.domain.type == "FISNA" & next.domain.type == "NACHT" & next.domain.start - current.domain.end +1 <= 1000]
  out.alldat[[7]] <- dt[previous.domain.type %in% c("B-BOX", "C3HC4") & current.domain.type == "fn3"]
  out.alldat[[8]] <- dt[(previous.domain.type %in% c("B-BOX", "C3HC4") & current.domain.type == "B30.2") | (previous.domain.type %in% c("B-BOX", "C3HC4") & current.domain.type == "fn3" & next.domain.type == "B30.2")]
  out.alldat[[9]] <- dt[current.domain.type == "NACHT" & next.domain.type == "B30.2"]
  out.alldat[[10]] <- dt[previous.domain.type == "NACHT" & current.domain.type == "FIIND" & next.domain.type == "CARD"]
  out.alldat[[11]] <- dt[previous.domain.type =="CARD" & (current.domain.type == "NACHT" | next.domain.type == "NACHT" & next.domain.start - previous.domain.end <= 99999)]
  out.alldat[[12]] <- dt[previous.domain.type =="PYRIN" & (current.domain.type == "NACHT" | next.domain.type == "NACHT" & next.domain.start - previous.domain.end <= 99999)]
  out.alldat[[13]] <- out.alldat[[12]][(current.domain.type == "NACHT" & current.domain.start %in% out.alldat[[9]]$current.domain.start) | (next.domain.type == "NACHT" & next.domain.start %in% out.alldat[[9]]$current.domain.start)]
  out.alldat[[14]] <- out.alldat[[11]][(current.domain.type == "NACHT" & current.domain.start %in% out.alldat[[9]]$current.domain.start) | (next.domain.type == "NACHT" & next.domain.start %in% out.alldat[[9]]$current.domain.start)]
  
  
  
  names(out.alldat) <- subsets
  
  
  out.counts <- c(assembly)
  dists.frames <- list()
  dists <- list()
  
  for (i in 1:14){
    out.counts <- c(out.counts, nrow(out.alldat[[i]]))
    tmpdat <- as.data.frame(table(out.alldat[[i]]$contig))
    if (length(tmpdat) != 2){
      tmpdat <- data.frame(Var1 = as.character("nope"), Freq = c("0"),stringsAsFactors=FALSE)
    }
    colnames(tmpdat) <- c("contig", subsets[i])
    dists[[i]] <- tmpdat
  }
  
  
  
  dists <- Reduce(function(x, y) merge (x, y, by = "contig", all=TRUE), dists)
  dists[is.na(dists)] <- 0
  dists <- dists[dists$contig != "nope" & order(as.character(dists$contig)),]
  
  cnt <- as.character(dists$contig)
  dists <- data.frame(lapply(dists, as.numeric, stringsAsFactors=FALSE))
  dists$contig <- cnt
  
  dists.nostrands <- dists
  dists.nostrands$contig <- gsub("_..$", "", dists.nostrands$contig)
  dists.nostrands <- aggregate(. ~ contig, dists.nostrands, sum)
  
  
  fwrite(dt, file=paste(assembly, ".dataset.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(dists, file = paste(assembly, ".contigs.strands.distrib.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(dists.nostrands, file = paste(assembly, ".contigs.distrib.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
  write.table(data.frame(out.counts), file = paste(assembly, ".counts.tsv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
  
  
  
  #follow it up in bash with the following: head -n1 GCA_900246225.contigs.distrib.tsv |cut -f2- | tr '\t' '\n' | sed -r '1s/^/\ncounted\n/' | paste - *counts.tsv | tr -d '\r' | sed '1d' | datamash transpose > counts.tsv; rm *.counts.tsv
}