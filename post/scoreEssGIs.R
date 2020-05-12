# Script to measure essential gene interactions (without qGI score)
library(plyr)
library(dplyr)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(gridExtra)
source("predGuides.R")

dataDir <- "/Users/catherineross/GIN"
dt <- format(Sys.Date(), "20%y%m%d")
outDir  <- sprintf("%s/res/%s_out_tiling_score_essGIs_predGuides", dataDir, dt)
if (!file.exists(outDir)) dir.create(outDir)

####################################################
# DATA INPUT
####################################################
## Essential genes
ess_tkoF <- sprintf("%s/data/pipeline_input/essentialGenes/essential_genes_6_12.txt", dataDir)
ess_tko <- read.delim(ess_tkoF, h=FALSE, as.is=TRUE)

## CRISPR scores (output from writeScoreTable.F)
val_scoreF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_val/table_tiling_guideScores.xlsx", dataDir)
val_score <- read.xlsx(val_scoreF)

#########################
## Tiling library data
#########################
# Dropout (foldchange) data
valF <- list.files(pattern="foldchange_mean.txt",
                   path=sprintf("%s/data/validation_tiling", dataDir),
                   full.names=TRUE)
val <- lapply(valF, function(x) read.delim(x, h=TRUE, as.is=TRUE))

# Filter to use only for effective guides (via predGuides.R)
guides <- predGuides(lib="tiling", dat=val_score, tko_guides=TRUE)

for (i in seq_along(val)) {
  val[[i]] <- val[[i]][which(val[[i]]$GENE_CLONE %in% guides),]
}

# Restucture data
## Query name info
str <- substr(basename(valF), 0, nchar(basename(valF))-29)
num <- substr(str, 4, 6)
query <- substr(str, 8, nchar(str))
name <- paste(query, num, sep="_")

## Gene and guide info
gene_clone_split <- strsplit(val[[1]]$GENE_CLONE, "_")
gene <- unlist(lapply(gene_clone_split, "[", 1))

## Put data together
fc <- do.call("cbind", lapply(val, "[", 4))
colnames(fc) <- name
fc <- as.matrix(fc)
rownames(fc) <- gene

## Subset for esssential genes
fc_ess <- fc[which(rownames(fc) %in% ess_tko$V1),]

####################################################
# PLOT LIBRARY WIDE LFC DISTRIBUTION
####################################################
# Plot guide LFC distribution per essential / non-essential gene

## NOTE indicate where tkov3 guides fall in distribution
plot_list <- list()
for (k in unique(rownames(fc))) {
  print(k)

  dat_plot <- fc[which(rownames(fc) == k),]
  dat_plot_melt <- melt(dat_plot)

  if (length(nrow(dat_plot)) == 0) {
    dat_plot_melt <- data.frame(Gene=k, Screen=rownames(dat_plot_melt), LFC=dat_plot_melt$value)
  }

  colnames(dat_plot_melt) <- c("Gene", "Screen", "LFC")

  if (k %in% ess_tko$V1) {
    title <- sprintf("%s guide dropout distribution (essential)", k)
  } else {
    title <- sprintf("%s guide dropout distribution (non-essential)", k)
  }

  if (length(nrow(dat_plot)) == 0) {
    subtitle <- "n=1 guide targetting gene"
  } else {
    subtitle <- sprintf("n=%s guides targetting gene", nrow(dat_plot))
  }

  plot_list[[k]] <-
    ggplot(dat_plot_melt, aes(x=LFC, fill=Screen)) +
      geom_density(alpha=0.2) +
      labs(x="Guide-level LFC", y="Density", title=title, subtitle=subtitle) +
      scale_fill_brewer(palette="Dark2") +
      theme_bw() +
      theme(text=element_text(family="sans", size=15))
}

plotF1 <- sprintf("%s/plot_tiling_distribution_guideLFC_perGene.pdf", outDir)
plot_grob <- marrangeGrob(plot_list, nrow=5, ncol=5)
ggsave(plotF1, plot_grob, width=40, height=35)

####################################################
# PLOT WT LFC DISTRIBUTION
####################################################
# To determine LFC threshold
pos_lim=-0.5
neg_lim=-1.5

fc_wt <- density(fc_ess[,"WT_001"])
plotF2 <- sprintf("%s/plot_tiling_density_wtLFC_essGenes.pdf", outDir)
pdf(plotF2, width=7, height=6)
plot(fc_wt, main="Distribution of guide-level essential gene LFC in WT background")
polygon(fc_wt, col="grey", border="black")
abline(v=c(pos_lim, neg_lim), col="red", lty=2, lwd=2)
dev.off()

####################################################
# MEASURE GIs
####################################################
# 1) Use the WT data to measure the mean LFC and variance for each individual
# guide in the tiling library
## NOTE only one wildtype screen in tiling data -- next step

# 2) For each essential gene, select subset of guides that have modest and
# stable LFC effects in WT screens
## Find guides that fall within specified [pos_lim, neg_lim] LFC range
#fc_sel <- fc_ess[which(fc_ess[,"WT_001"] <= pos_lim & fc_ess[,"WT_001"] >= neg_lim),]

# Get number of guides per gene and write out
#n_guide <- as.data.frame(table(rownames(fc_sel)))
#n_guide <- n_guide[order(n_guide$Freq),]
#tmpF <- sprintf("%s/table_numGuideGene2.txt", outDir)
#write.table(n_guide, tmpF, col=TRUE, row=FALSE, quote=F, sep="\t")

# Get genes that are excluded after filtering, if any
#genes_org <- unique(rownames(fc_ess))
#genes_sel <- unique(rownames(fc_sel))
#matches <- which(genes_org %in% genes_sel)
#genes_org[-matches]

# 3) Measure the diff LFC (Query LFC - WT LFC) for the selected guides,
# take the mean and SD of DIFF LFC to call GIs for each query
# NOTE using predGuides set, not using LFC to filter guides
fc_sel <- fc_ess
fc_diff <- matrix(nrow=nrow(fc_sel), ncol=5, dimnames=list(rownames(fc_sel), name[2:6]))

for (i in 2:ncol(fc_sel)) {
  diff <- fc_sel[,i] - fc_sel[,1]
  fc_diff[,i-1] <- diff
}

# Sumarize diff LFC by gene
gene_unique <- unique(rownames(fc_diff))
fc_diff_mean <- matrix(nrow=length(gene_unique), ncol=ncol(fc_diff), dimnames=list(gene_unique, colnames(fc_diff)))
fc_diff_sd <- matrix(nrow=length(gene_unique), ncol=ncol(fc_diff), dimnames=list(gene_unique, colnames(fc_diff)))

for (k in seq_along(gene_unique)) {
  print(k)
  gene_k <- gene_unique[k]
  gene_fc <- fc_diff[which(rownames(fc_diff) == gene_k),]
  # Calculate mean per gene
  if (length(nrow(gene_fc)) == 0) {
    fc_diff_mean[k,] <- gene_fc
    fc_diff_sd[k,] <- 0
  } else {
    gene_fc_mean <- colMeans(gene_fc, na.rm=TRUE)
    fc_diff_mean[k,] <- gene_fc_mean
    # Calculate sd per gene
    gene_fc_sd <- apply(gene_fc, 2, sd, na.rm=TRUE)
    fc_diff_sd[k,] <- gene_fc_sd
  }
}

# Write out table of diff LFC mean
outF <- sprintf("%s/table_tiling_scores_LFC_essGenes.txt", outDir)
write.table(fc_diff_mean, outF, col=TRUE, row=TRUE, quote=FALSE, sep="\t")

# Organize by query and write out
fc_out <- list()

for (x in 1:ncol(fc_diff_mean)) {
  name_x <- colnames(fc_diff_mean[,x,drop=FALSE])
  dat_mean <- fc_diff_mean[,x]
  dat_sd <- fc_diff_sd[,x]
  dat_comb <- cbind(DIFF_LFC_MEAN=dat_mean, DIFF_LFC_SD=dat_sd)
  dat_comb <- dat_comb[order(dat_comb[,1]),]
  fc_out[[name_x]] <- dat_comb
}

outF2 <- sprintf("%s/table_tiling_scores_LFC_essGenes.xlsx", outDir)
write.xlsx(fc_out, outF2, row.names=TRUE)
