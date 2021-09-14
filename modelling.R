# If using corPagel for modelling below, keep in mind the following.
# Lambda = 0 -> phylogeny has no effect. Same as not including a tree into the model.
# Lambda = 1 -> phylogeny has the maximal effect and has a Brownian distribution. Same as using corBrownian

setwd("E:/Dropbox/TRIMpaper_2021/with_RING-NLR")
Sys.setenv(ENTREZ_KEY="") # insert your key here, needed to access NCBI databases

library(taxize)
library(ape)
library(geiger)
library(nlme)
library(ggplot2)
library(sf)
library(rnaturalearth)



dfdat <- read.table("E:/Dropbox/TRIMpaper_2021/onlyfish_values_v11.tsv", sep = "\t", header = T)
domtypes <- colnames(dfdat)[8:21]

# dfdat <- dfdat[dfdat$full_chroms == 1 & dfdat$Long_reads == 1,]
dfdat <- dfdat[dfdat$full_chroms == 1 & dfdat$Long_reads == 1,]

# create a new categorical variable for water/habitat: -1 = freshwater, 0 = either/both, 1 = marine
dfdat$water <- dfdat$marine - dfdat$freshwater

# convert latitude values to distances from the equator
dfdat$lat_c <- abs(dfdat$lat_c)

# plot coordinates on a world map

coordf <- st_as_sf(dfdat[c("lat_c", "water")], coords = c("water", "lat_c"), crs = "WGS84")
world <- ne_countries(scale = "medium", returnclass = "sf")
pdf("worldmap.pdf", width = 11.6, height = 7.7)
ggplot(data=world) + geom_sf(fill = "antiquewhite1") + geom_sf(data=coordf) + coord_sf(xlim=c(-180, 180), ylim=c(-90, 90), expand = FALSE) + geom_vline(xintercept = 0) + geom_hline(yintercept=0) + xlab("longitude") + ylab("latitude")
dev.off()


# make a tree using Species names
Species <- dfdat$Species

# get NCBI unique identifier (UID) for each Species:
uids <- lapply(as.list(Species), function(x) get_uid(x, messages = FALSE)[1])

# append UID to data because why not:
dfdat$uid <- unlist(uids)

# make Species names rownames for the model
rownames(dfdat) <- Species

# create a tree:
taxize_class <- classification(uids, db = "ncbi")
taxize_tree <- class2tree(taxize_class, check = TRUE)
fishtree <- taxize_tree$phylo # this is the class "phylo" object that goes in the model
plot(taxize_tree)


# Check to make sure names match up. If some don't fix them.
name.check(fishtree, dfdat)

# calculate and plot the likelihood surface of the model in response to different lambda values and different domain categories

fishplots <- list()
lambdadat <- list()
lambda <- seq(0, 1, 0.001)

pagel <- function(domtype, lambda){
  gls(as.formula(paste("scale(", domtype, ") ~ scale(Genome_len) + scale(lat_c) + scale(water)")),
             correlation = corPagel(lambda, fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude)
}

for (domtype in domtypes){
  lik <- sapply(lambda, function(lambda) logLik(pagel(domtype, lambda)))
  lambdalik <- data.frame(cbind(lik,lambda))
  lambdadat[[domtype]] <- lambdalik
  fishplots[[domtype]] <- ggplot(lambdalik, aes(x = lambda, y = lik)) + geom_line() +
    theme_bw() + ggtitle(paste0(domtype, " in all fish")) + xlab("lambda") + ylab("log likelihood")
}

pdf("lambdaplots.pdf")
fishplots
dev.off()


# define a function to estimate lambda automatically

pagel_estimate_lambda <- function(domtype){
  gls(as.formula(paste("scale(", domtype, ") ~ scale(Genome_len) + scale(lat_c) + scale(water)")),
      correlation = corPagel(1, fishtree, form = ~Species), data = dfdat, method = "ML", na.action=na.exclude)$modelStruct
}




# NB! NOT AUTOMATED! Next run pagel_estimate_lambda for each domain type to get the lambda values resulting in best model fits.
# if lambda estimation fails, calculate the associated likelihood surface for the particular domain type, choose the likeliest value, and do a sanity check by checking the plot
# note that lambda values less than 0 or more than 1 are artifacts, CHECK THE PLOT, it is probably just 0 or 1
# example: llik <- sapply(lambda, function(lambda) logLik(pagel("B.Box", lambda)))
# the following values are from either lambda estimation OR checking for one producing likeliest models, as described above

bestfits <- data.frame(domtypes, c(0.9004214, 0.892, 0.9454092, 0.9591276, 0.7261048, 0.794984, 0.8071393, 0.9250062, 0.7983087, 1, 0.9951094, 0.6552412, 0.625054, 0.8487437))
colnames(bestfits) <- c("domtype", "lambda")


# run the models and dump all summaries into a file. This will then have to be processed into a proper table before proceeding.
sink('model_summaries.tsv')
summary(gls(scale(B30.2.PRY.SPRY) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "B30.2.PRY.SPRY", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(B.Box) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "B.Box", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(TRIM) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "TRIM", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(TRIM.B30.2) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "TRIM.B30.2", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(TRIM.fn3) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "TRIM.fn3", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(NLR) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "NLR", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(FISNA.NACHT) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "FISNA.NACHT", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(RING.NLR) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "RING.NLR", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(NLR.B30.2) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "NLR.B30.2", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(PYD.NLR) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "PYD.NLR", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(PYD.NLR.B30.2) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "PYD.NLR.B30.2", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(CARD.NLR) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "CARD.NLR", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(CARD.NLR.B30.2) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "CARD.NLR.B30.2", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
summary(gls(scale(NLR.FIIND.CARD) ~ scale(Genome_len) + scale(lat_c) + scale(water), correlation = corPagel(bestfits[bestfits$domtype == "NLR.FIIND.CARD", 2], fishtree, form = ~Species, fixed = TRUE), data = dfdat, method = "ML", na.action=na.exclude))
sink()


# NB! NOT AUTOMATED! Convert the file model_summaries.tsv into coefficient tables, etc, with the script "tabulate_summaries.sh" from bash


# The below is just a reminder of what should be in the coefficients file before loading them back to R. The others are not used in this version of the script.
# colnames(df_coef) <- c("domains", "param", "value", "stderr", "tval", "pval")
# colnames(df_intercept) <- c("domains", "value", "stderr", "tval", "pval")
# colnames(df_likelihood) <- c("domains", "AIC", "BIC", "loglik")
# colnames(df_residuals) <- c("domains", "Min", "Q1", "Med", "Q3", "Max", "stderr")


# read in the output from tabulate_summaries.sh

df_coef <- read.table("model_coefficients.tsv", stringsAsFactors = FALSE, header = FALSE)

# rename columns, exclude domain types for which plotting is not needed

colnames(df_coef) <- c("domains", "param", "value", "stderr", "tval", "pval")
df_coef <- df_coef[! df_coef$domains == "TRIM.B.Box",]


# calculate 95% confidence intervals

df_coef$stderr_low <- df_coef$value - 1.96*df_coef$stderr
df_coef$stderr_high <- df_coef$value + 1.96*df_coef$stderr


# transform the column used for labelling into a factor to stop ggplot from reordering things

df_coef$domains <- factor(df_coef$domains, levels=rev(unique(df_coef$domains)))


#plot the data

library("gridExtra")

plot.coefficient.densities <- function(dataset, datalabel){
  prettyplot <- function(x, y){
    ggplot(x,aes(x=domains, y=value)) +
      geom_rect(xmin = 9.5, xmax = 13.5, ymin = -Inf, ymax = +Inf, fill = "Gainsboro", alpha = 0.5) +
      geom_hline(yintercept = 0) + 
      geom_point(shape = 3) + 
      geom_linerange(aes(x = domains, ymin = stderr_low, ymax = stderr_high)) + 
      ylim(-1, 2) +
      theme_bw() +
      theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 1, hjust=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_text(size = 8, aes(y = stderr_high, label=ifelse(pval < 0.05, "*", ""))) +
      geom_text(size = 8, aes(y = stderr_high+0.15, label=ifelse(pval < 0.01, "*", ""))) +
      geom_text(size = 8, aes(y = stderr_high+0.3, label=ifelse(pval < 0.001, "*", ""))) +
      ylab("Effect size") +
      ggtitle(y) +
      coord_flip()
  }
  
  
  gensize <- dataset[dataset$param == "Genome_len",]
  gensize_sig <- gensize[gensize$pval < 0.05,]
  latc <- dataset[dataset$param == "lat_c",]
  latc_sig <- latc[latc$pval < 0.05,]
  wat <- dataset[dataset$param == "water",]
  wat_sig <- wat[wat$pval < 0.05,]
  
  Aplot <- prettyplot(gensize, "Genome size")
  Bplot <- prettyplot(latc, "Distance from equator")
  Cplot <- prettyplot(wat, "Habitat")
  
  
  grid.arrange(Aplot, Bplot, Cplot, ncol = 3, nrow = 1)
}


pdf("coefficient_distribs_onlychroms_v3.pdf", width = 10, height = 5)
plot.coefficient.densities(df_coef, "predicted gene counts")
dev.off()



# create additional plots to explore some of the domain types

df1 <- dfdat[,c(6,11,5,26)]
df2 <- dfdat[,c(6,16,5,28)]
df1$domtype <- "TRIM.B30.2"
df2$domtype <- "NLR.B30.2"
colnames(df1) <- c("latitude", "val", "genome_len", "clustering", "domtype" )
colnames(df2) <- colnames(df1)
df1$scval <- df1$val/max(df1$val)
df2$scval <- df2$val/max(df2$val)
df <- rbind(df1, df2)

plot1 <- ggplot(df, aes(domtype, val)) +
  geom_violin(aes(fill = domtype)) +
  scale_x_discrete(labels= c("NLR-B30.2", "TRIM-B30.2")) +
  scale_fill_manual(values = c("Gainsboro", "DarkGray")) +
  theme_bw() + theme(axis.title.x=element_blank(), legend.position = "none") +
  ylab("Counts")

plot2 <- ggplot(df, aes(domtype, clustering)) +
  geom_violin(aes(fill = domtype)) +
  scale_x_discrete(labels= c("NLR-B30.2", "TRIM-B30.2")) +
  scale_fill_manual(values = c("Gainsboro", "DarkGray")) +
  theme_bw() + theme(axis.title.x=element_blank(), legend.position = "none") +
  ylab("Min # of scaffolds combined for 25%")

plot3 <- ggplot(dfdat, aes(NLR.B30.2, TRIM.B30.2)) +
  geom_point() +
  theme_bw() + theme(legend.title = element_blank()) +
  xlab("NLR-B30.2 counts") + ylab("TRIM-B30.2 counts")

plot4 <- ggplot(dfdat, aes(L25_NLR.B30.2, L25_TRIM.B30.2)) +
  geom_point() +
  xlim(0,12) +
  theme_bw() + theme(legend.title = element_blank()) +
  xlab("Min # of scaffolds combined for 25% of TRIM-B30.2") + ylab("Min # of scaffolds combined for 25% of NLR-B30.2")

plot5 <- ggplot(df, aes(scval, genome_len)) +
  geom_point(aes(color = domtype)) +
  scale_shape_manual(values=c(16, 3)) + scale_color_manual(values = c("black", "red")) +
  labs(x="Domain counts (scaled to 1)", y = "Genome size") +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom") +
  scale_y_continuous(labels= c("1 Gb", "2 Gb", "3 Gb"), breaks = c(1000000000, 2000000000, 3000000000))

plot6 <- ggplot(df, aes(scval, latitude)) +
  geom_point(aes(color = domtype)) +
  scale_shape_manual(values=c(16, 3)) + scale_color_manual(values = c("black", "red")) +
  labs(x="Domain counts (scaled to 1)", y = "Distance from the equator (degrees)") +
  theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom")

plot7 <- ggplot(dfdat, aes(as.character(water), NLR.B30.2)) +
  geom_violin(aes(fill = water)) +
  scale_x_discrete(labels= c("freshwater", "both", "marine")) +
  theme_bw() + theme(axis.title.x=element_blank(), legend.position = "none")

plot8 <- ggplot(dfdat, aes(as.character(water), TRIM.B30.2)) +
  geom_violin(aes(fill = water)) +
  scale_x_discrete(labels= c("freshwater", "both", "marine")) +
  theme_bw() + theme(axis.title.x=element_blank(), legend.position = "none")

grid.arrange(plot1, plot2, plot5, plot3, plot4, plot6, plot7, plot8, ncol = 3)
