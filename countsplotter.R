setwd("/Dropbox/TRIMpaper_2021/with_RING-NLR")

library(plyr)
library(reshape)
library(ggplot2)
library(scales)


df <- read.table("/Dropbox/TRIMpaper_2021/targetfish_v11_for_plotting.tsv", sep = "\t", header = T)
df$Species <- paste(df$Species, df$order, sep  = " - ")
species_df <- df$Species
df$Genome_length <- round(df$Genome_length/1000000000,1)
df$N50 <- round(df$N50/1000000,2)

# do the following only once, the goal is to get a tree, then use the tree order to re-sort the table targetfish_v10.tsv
# uids_df <- lapply(as.list(species_df), function(x) get_uid(x, messages = FALSE)[1])
# taxize_class_df <- classification(uids_df, db = "ncbi")
# taxize_tree_df <- class2tree(taxize_class_df, check = TRUE)
# pdf("allspecies.pdf", height = 30)
# plot(taxize_tree_df)
# dev.off()

df2 <- df[,c(1,3:6,9:12,14:20,49,50)]
df2 <- melt(df2)
df2<- ddply(df2, .(variable), transform, rescale = rescale(value))
df2$Species <- factor(df2$Species, levels = rev(species_df))

pdf("heatmap.pdf", height = 30, width = 10)
ggplot(df2, aes(variable, Species)) + 
  geom_tile(aes(fill = rescale), colour = "white") +
  geom_text(aes(label=value)) + 
  scale_fill_gradient(low = "white", high = "DarkSlateBlue", na.value = "black") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none") +
  scale_x_discrete(position = "top")
dev.off()