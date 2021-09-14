### finds the minimum number of contigs containing a percentage of the domain of interest

setwd(system("pwd", intern = T) )


args <- commandArgs(TRUE)

infile.name <- args[1]
n.percent <- as.numeric(args[2])/100

df <- read.table(args[1], header = T)
domvec <- infile.name



for (domstruct in 2:15) {
  a <- sort(as.numeric(df[,domstruct]), decreasing = TRUE)
  
  counter <- 0
  ndom <- 0
  for (i in a) {
    if (ndom < sum(a) * n.percent ){
      counter <- counter + 1
      ndom <- ndom + i
    }
  }
  domvec[domstruct] <- counter
}

cat(domvec)


# example usage: rm *n25.tsv *n50.tsv; date; for i in *contigs.distrib.tsv; do Rscript contig-counter.R $i 25 >> contigs.distribs.n25.tsv; Rscript contig-counter.R $i 50 >> contigs.distribs.n50.tsv; done; date
