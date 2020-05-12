# Script to select essential genes for query construction
## 1. Rank essential genes based on their array side negative interaction degree
##    (standard cutoff).
## 2. Focus on essential genes that have high negative degree and also show
##    negative interactions with 2/4 replicates of either PDCD5, FASN or FANCG
## 3. #2 is a strict filter so if we don’t have any essentials that interact
##    reproducibly with at least 2 of those 4 queries then we’ll expand the
##    analysis to include all queries screened thus far
library(openxlsx)
library(plyr)
library(dplyr)
library(OneR)
library(fgsea)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)

dataDir <- "/Users/catherineross/GIN"
outDir  <- sprintf("%s/bin/R/guideAnalysis/output_data/out_TKOv3/plot_hypomorph", dataDir)

####################################################
# DATA INPUT
####################################################
## Experimental info
#########################
## Essential genes
### TKOv3
ess_tkoF <- sprintf("%s/data/pipeline_input/essentialGenes/essential_genes_6_12.txt", dataDir)
ess_tko <- read.delim(ess_tkoF, h=FALSE, as.is=TRUE)

### Tiling
ess_valF <- sprintf("%s/data/validation_tiling/tiling_essentials.txt", dataDir)
ess_val <- read.delim(ess_valF, h=FALSE, as.is=TRUE)

## Pathway annotations
annoF <- sprintf("%s/anno/Human_GO_bp_no_GO_iea_October_01_2018_symbol.gmt", dataDir)
anno <- gmtPathways(annoF)
anno_len <- lapply(anno, length)
anno2 <- anno[-which(anno_len >= 500)] # remove pathways with >500 genes

#########################
## TKOv3 qGI data
#########################
## NOTE using most recent qGI data
inDir <- sprintf("%s/bin/R/huGI/20191015_automated/output_data/20191023_replicateMeanPairs_mergeFC", dataDir)
pisF <- list.files(pattern="PiScores", path=inDir, full.names=TRUE)
pis  <- read.delim(pisF, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
fdrF <- list.files(pattern="FDR", path=inDir, full.names=TRUE)
fdr  <- read.delim(fdrF, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)
# Remove data for 'comb' screens
pis <- pis[,-grep("comb", names(pis))]
fdr <- fdr[,-grep("comb", names(fdr))]

#########################
## Tiling library data
#########################
# Dropout (foldchange) data
valF <- list.files(pattern="foldchange_mean.txt",
                   path=sprintf("%s/data/validation_tiling", dataDir),
                   full.names=TRUE)
val <- lapply(valF, function(x) read.delim(x, h=TRUE, as.is=TRUE))

# Merge query foldchange data together (expected all behave similarly)
val_all <- matrix(nrow=nrow(val[[1]]), ncol=length(val))
for (i in 1:ncol(val_all)) {
  val_all[,i] <- val[[i]]$T19_Minimal
}
val_mean <- as.data.frame(rowMeans(val_all, na.rm=TRUE))
val_mean <- cbind(guide=val[[1]]$GENE_CLONE, gene=val[[1]]$GENE, val_mean)
colnames(val_mean)[3] <- "logFC"
val_mean$guide <- gsub(".*\\_", "", val_mean$guide) # get only guide sequence

# Get mean logFC per gene
val_gene_mean <- aggregate(val_mean$logFC, list(val_mean$gene), mean)
val_gene_mean <- na.omit(val_gene_mean)
colnames(val_gene_mean) <- c("gene", "mean_logFC")

#########################
## GIN wildtype data
#########################
# Wt dropout (foldchange) data
wtF <- list.files(pattern="GIN.*.xlsx", path=sprintf("%s/data/wildtypes", dataDir), full.names=TRUE)
wt_name <- paste("wt", substr(basename(wtF), 4, 6), sep="_")
wt <- lapply(wtF, function(x) read.xlsx(x, sheet=2))

for (i in seq_along(wt)) {
  tmp <- wt[[i]]
  colnames(tmp)[4:ncol(tmp)] <- paste(wt_name[i], colnames(tmp)[4:ncol(tmp)], sep="_")
  wt[[i]] <- tmp
}

# Join data and keep relevant columns
to_rm <- paste(c("Rich", "pyruvate", "PTC"), collapse="|")
tps <- paste(c("T18", "T19", "T20"), collapse="|")
wt_all <- join_all(wt, by=c("Guide.ID", "GENE_CLONE", "GENE"))
wt_dat <- wt_all[,-grep(to_rm, colnames(wt_all))]
wt_dat <- wt_dat[,grep(tps, colnames(wt_dat))]
colnames(wt_dat) <- substr(colnames(wt_dat), 0, nchar(colnames(wt_dat))-12)

# Get mean wt log2FC
wt_mean <- as.data.frame(rowMeans(wt_dat))
wt_mean <- cbind(guide=wt_all$GENE_CLONE, gene=wt_all$GENE, wt_mean)
colnames(wt_mean)[3] <- "logFC"
wt_mean$guide <- gsub(".*\\_", "", wt_mean$guide) # get only guide sequence

# Get mean logFC per gene
wt_gene_mean <- aggregate(wt_mean$logFC, list(wt_mean$gene), mean, na.rm=TRUE)
wt_gene_mean <- na.omit(wt_gene_mean)
colnames(wt_gene_mean) <- c("gene", "mean_logFC")

####################################################
# DATA ANALYSIS
####################################################
# Filter for essential genes
#ess_set = ess_val
#pis_ess <- pis[which(rownames(pis) %in% ess_set[,1]),]
#fdr_ess <- fdr[which(rownames(fdr) %in% ess_set[,1]),]

# Get mean logFC for tiling essentails
ess_val_mean <- val_gene_mean[which(val_gene_mean$gene %in% rownames(pis_ess)),]

# Rank by neg GI degree
ranks <- data.frame(gene=wt_gene_mean$gene, meanLogFC_wtTKO=wt_gene_mean$mean_logFC,
                    negDen=NA, repGI=NA, repNum_PDCD5=NA, repNum_FASN=NA, repNum_FANCG=NA)
ranks$gene <- as.character(ranks$gene)

# Focus on high neg degree genes that interact with replicates
reps <- paste(c("PDCD5", "FASN", "FANCG"), collapse="|")

for (i in 1:nrow(ranks)) {
  print(ranks[i,]$gene)
  gene <- ranks[i,]$gene
  dat <- pis[gene, which(pis[gene,] <= -0.5 & fdr[gene,] <= 0.5)]
  queries <- names(dat)
  queries_rep <- queries[grep(reps, queries)]

  ranks[i,]$negDen <- round((length(dat)/ncol(pis))*100, 2)
  ranks[i,]$repGI <- paste(sort(queries_rep), collapse=", ")
  ranks[i,]$repNum_PDCD5 <- length(grep("PDCD5", queries_rep))
  ranks[i,]$repNum_FASN <- length(grep("FASN", queries_rep))
  ranks[i,]$repNum_FANCG <- length(grep("FANCG", queries_rep))
}

# Sort by decreasing order
ranks <- ranks[order(ranks[,2], ranks[,3], decreasing=FALSE),]

# Filter out genes with sig neg GI in 0 / 1 replicate query
no_sig <- which(ranks[,"repNum_PDCD5"] <= 1 & ranks[,"repNum_FASN"] <= 1 & ranks[,"repNum_FANCG"] <= 1)
ranks_filt <- ranks[-no_sig,]

# Also filter out essential genes that have been screened as a query
#query <- substr(names(pis_ess), 0, nchar(names(pis_ess))-4)
#query_gene <- which(query %in% ranks_filt$gene)
#if (!length(query_gene)) {
#  ranks_filt2 <- ranks_filt
#} else {
#  ranks_filt2 <- ranks_filt[-query_gene,]
#}

# Put all data into list and write out
ranks_list <- list()
ranks_list[["filtered"]] <- ranks_filt
ranks_list[["all"]] <- ranks

# Write out excel sheet
outF <- sprintf("%s/table_TKOv3_minWTlogFC_negDen_repGI_GIN20191015.xlsx", outDir)
write.xlsx(ranks_list, outF, row.names=FALSE)

####################################################
# ADDITIONAL ANALYSIS
####################################################
# Read in excel sheet with Julianne annotation of Horizon KO genes
dat_koF <- sprintf("%s/table_TKOv3_essGenes_negDeg_repGI_GIN20191015_JL.xlsx", outDir)
dat_ko <- read.xlsx(dat_koF)
dat_ko <- dat_ko[,c(1,2)]

# Filter from table generated for tiling essential genes
dat_ko_val <- merge(dat_ko, ranks_filt2, by="gene")
dat_ko_val <- dat_ko_val[order(dat_ko_val$meanLogFC_val),]

# Select those without KO lines
dat_noKO_val <- dat_ko_val[which(dat_ko_val[,2] == "no"),]

# Add the tiling guides + LFC per gene
dat_noKO_val_guides <- merge(val_mean, dat_noKO_val, by="gene")
dat_noKO_val_guides <- dat_noKO_val_guides[order(dat_noKO_val_guides$meanLogFC_val, dat_noKO_val_guides$logFC),]

# Write out
outF2 <- sprintf("%s/table_tiling_essGenes_negDen_repGI_noKO_GIN20191015.xlsx", outDir) # used ess_val set
write.xlsx(dat_noKO_val_guides, outF2, row.names=FALSE)

####################################################
# Updated MC essential gene candidate list with proper degree / density values
glistF <- sprintf("%s/essential gene candidates.xlsx", outDir)
glist  <- read.xlsx(glistF)
degF <- sprintf("%s/res/20191115_out_networkGIdensity_20191023_replicateMeanPairs_mergeFC/networkDensity_whole_standard_singleStats_neg.txt", dataDir)
deg  <- read.delim(degF, h=TRUE, as.is=TRUE)
# Fill in data for genes not in table (those with degree=0)
deg <- na.omit(deg)
toFill <- rownames(pis)[-which(rownames(pis) %in% deg$Var1)]
toFill <- data.frame(Var1=toFill, Freq=0, density=0)
deg <- rbind(deg, toFill)

for (i in 1:nrow(glist)) {
  print(i)
  gene <- glist[i,]$gene
  deg_replace <- deg[which(deg$Var1 == gene),]$Freq
  den_replace <- deg[which(deg$Var1 == gene),]$density

  if (!length(deg_replace)) {
    glist[i,]$neg.degree <- 0
  } else {
    glist[i,]$neg.degree <- deg_replace
  }

  if (!length(den_replace)) {
    glist[i,]$neg.density <- 0
  } else {
    glist[i,]$neg.density <- round(den_replace, 2)
  }
}

# Incorporate guide behaviour (TKOv3 guide LFC) data
glist_tko_guides <- merge(wt_mean, glist, by="gene")
glist_tko_guides <- glist_tko_guides[order(glist_tko_guides$meanLogFC_wtTKO, glist_tko_guides$logFC),]

outF3 <- sprintf("%s/table_essCandidates_guideLFC.xlsx", outDir)
write.xlsx(glist_tko_guides, outF3)

####################################################
# GO BP ANALYSIS
####################################################
df <- read.xlsx(outF3)

# List relevant bioprocesses per gene
df$GOBP <- NA
for (i in unique(df$gene)) {
  print(i)
  gene = sprintf("^%s$", i)
  gene_anno <- lapply(anno2, function(x) grep(gene, x))
  gene_anno <- names(which(gene_anno != 0))
  gene_anno <- gsub("\\%.*", "", gene_anno)
  df[which(df$gene == i),]$GOBP <- paste(gene_anno, collapse="; ")
}

outF4 <- sprintf("%s/table_essCandidates_guideLFC_GOBP.xlsx", outDir)
write.xlsx(df, outF4)

####################################################
# INTEGRATION OF CRISPR TOOLS
####################################################
df <- read.xlsx(outF4) # essential gene candidates

# Azimuth predictions on TKOv3 guides
azimuthF <- list.files(pattern="azimuth", path=outDir, full.names=TRUE)
azimuth <- read.delim(azimuthF, h=TRUE, stringsAsFactors=FALSE)
azimuth <- azimuth[,-3]
colnames(azimuth) <- c("gene", "guide", "azimuth_guideEfficiency_pred")

# Merge with table
df_pred <- join(df, azimuth)
outF5 <- sprintf("%s/table_essCandidates_guideLFC_GOBP_predictions.xlsx", outDir)
write.xlsx(df_pred, outF5)

####################################################
# DATA PREP FOR PLOTS
####################################################
# Integrate with guide score data (crispro, azimuth, indelphi)
tkoF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_TKOv3/plot_hypomorph/table_TKOv3_guideScores.xlsx", dataDir)
valF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_val/plot_hypomorph/table_tiling_guideScores.xlsx", dataDir)
tko <- read.xlsx(tkoF)
val <- read.xlsx(valF)

# Dataset to plot (choose one)
## All guide score data
dat_merge <- rbind(tko, val)
dat_merge$gene_fraction_bin <- bin(dat_merge$gene_fraction, nbins=10, method="length")
colnames(dat_merge)[2] <- "gene"

## Latest essential gene candidate table
df <- read.xlsx(outF5)

## Negative density table for all TKOv3 library genes
df <- deg
colnames(df) <- c("gene", "neg.degree", "neg.density")
## Filter for essentials
df <- df[which(df$gene %in% ess_tko$V1),]
# Incorporate guide behaviour (TKOv3 guide LFC) data
df <- join(df, wt_gene_mean, by="gene")

####################################################
# PLOTS
####################################################
## 1) GI degree vs. gene level dropout
outF <- sprintf("%s/point_TKOv3_ess_GIdeg_geneLFC.pdf", outDir)
pdf(outF, height=6, width=9)
p1 <- ggplot(df, aes(x=neg.degree, y=mean_logFC)) +   # meanLogFC_wtTKO
        geom_point() +
        geom_smooth() +
        theme_bw() +
        labs(x="Negative GI degree", y="Gene-level LFC",
             title="Relationship between essential gene negative GI degree and WT LFC",
             subtitle=sprintf("GI degree vs. gene-level dropout (N=%s genes)", nrow(df))) +
        theme(text=element_text(family="sans", size=15))
print(p1)
dev.off()

## 2) GI degree vs. guide level dropout
outF2 <- sprintf("%s/point_TKOv3_essCandidates_GIdeg_guideLFC.pdf", outDir)
pdf(outF2, height=6, width=9)
p2 <- ggplot(df, aes(x=neg.degree, y=logFC)) +
        geom_point() +
        geom_smooth() +
        theme_bw() +
        labs(x="Negative GI degree", y="Guide-level LFC",
             title="Relationship between essential gene negative GI degree and WT LFC",
             subtitle="GI degree vs. guide-level dropout") +
        theme(text=element_text(family="sans", size=15))
print(p2)
dev.off()

## 3) mean LFC vs. gene position (exon / gene fraction)
plot_list <- list()
for (k in unique(df$gene)) {
  print(k)
  dat_plot <- filter(dat_merge, gene == k)
  stat <- "gene_fraction"
  title <- k

  # Get gene-level mean LFC in each dataset
  means = ddply(dat_plot, ~dataset, summarise, mean=mean(logFC))
  dat_plot$geneLFC <- NA
  dat_plot[which(dat_plot$dataset == "TKOv3"),]$geneLFC <- means[which(means$dataset == "TKOv3"),]$mean
  dat_plot[which(dat_plot$dataset == "Tiling"),]$geneLFC <- means[which(means$dataset == "Tiling"),]$mean

  plot_list[[k]] <-
    ggplot(dat_plot, aes(x=factor(gene_fraction_bin), y=logFC, colour=guide_library, shape=guide_target)) +
        facet_wrap(.~dataset) +
        geom_point(size=2, alpha=0.5) +
        geom_errorbar(aes(ymin=logFC-sd, ymax=logFC+sd), width=0.2, alpha=0.5) +
        scale_colour_manual(values=c("black", "red")) +
        geom_hline(aes(yintercept=geneLFC), linetype="dashed", colour="blue") +
        labs(x=stat, y="Mean guide-level LFC\n(min wildtype TKOv3)", title=title,
             subtitle=sprintf("Guide LFC distribution vs %s (y=gene LFC)", stat)) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
}

outF3 <- sprintf("%s/point_TKOv3_tiling_essCandidates_guideLFC_vs_%s_perGene.pdf", outDir, stat)
plot_grob <- marrangeGrob(plot_list, nrow=5, ncol=3)
ggsave(outF3, plot_grob, width=40, height=35)

## 4) guide LFC consistency across libraries
tko_fc <- tko[,c("guide", "gene_name", "logFC")]
colnames(tko_fc)[3] <- paste(colnames(tko_fc)[3], "TKO", sep="_")
val_fc <- val[,c("guide", "gene_name", "logFC")]
colnames(val_fc)[3] <- paste(colnames(val_fc)[3], "tiling", sep="_")
# Combine logFC data for guides in both TKOv3 and tiling and get difference
fc_comb <- join(tko_fc, val_fc)
fc_comb <- na.omit(fc_comb)
# Subset for essential gene candidates
fc_comb <- fc_comb[which(fc_comb$gene_name %in% df$gene),]
fc_comb$fc_diff <- abs(fc_comb$logFC_TKO - fc_comb$logFC_tiling)
fc_comb <- fc_comb[order(fc_comb$fc_diff, decreasing=TRUE),]
fc_comb$int <- 1:nrow(fc_comb)
# Label guides >1 foldchange difference
fc_comb$label <- NA
toLabel <- which(fc_comb$fc_diff >= 1)
fc_comb[toLabel,]$label <- paste(fc_comb[toLabel,]$gene_name, fc_comb[toLabel,]$guide, sep="_")

## Plot difference in LFC measurements across libraries
outF4 <- sprintf("%s/point_TKOv3_tiling_essCandidates_LFC_consistency.pdf", outDir)
pdf(outF4, width=9, height=5)
p4 <- ggplot(fc_comb, aes(x=int, y=fc_diff, label=label)) +
        geom_point(size=1.5, alpha=0.5) +
        labs(x="Guide (in both libraries)", y="Guide-level LFC difference\n(TKOv3 vs. tiling)",
             title="LFC consistency for essential gene candidate guides",
             subtitle=sprintf("N=%s genes and %s guides (y=mean difference)", length(unique(fc_comb$gene)), nrow(fc_comb))) +
        geom_hline(aes(yintercept=mean(fc_diff)), linetype="dashed", colour="blue") +
        geom_text_repel(parse=TRUE, na.rm=TRUE, segment.size=0.2, segment.color="grey50",
                        nudge_x=500-subset(fc_comb, label != "NA")$int,
                        size=3, direction="y", hjust=0) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
print(p4)
dev.off()

## Correlation between logFC and gRNA efficiency (Azimuth prediction)
df_plot <- dat_merge # all TKOv3/tiling guides
df_plot <- df        # essential gene candidates

x = df_plot$azimuth_guideEfficiency_pred
y = df_plot$logFC

lm_eqn <- function(x, y) {
    m <- lm(y ~ x);
    r2 <- format(summary(m)$r.squared, digits = 3)
    eq <- sprintf("R2 = %s", r2)
    return(eq)
}

outF5 <- sprintf("%s/point_TKOv3_tiling_all_guideLFC_vs_gRNA_efficiency.pdf", outDir)
pdf(outF5, width=10, height=5)
p5 <- ggplot(df_plot, aes(x=x, y=y)) +
        facet_wrap(.~dataset, scales="free_x") +
        geom_point(alpha=0.3) +
        geom_smooth(method="lm") +
        labs(x="Predicted gRNA efficiency", y="Guide-level LFC",
             title="Guide-level dropout vs. gRNA efficiency predictions") +
        #geom_text(x=0.25, y=2, label=lm_eqn(x, y)) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15))
print(p5)
dev.off()

### Combine all plots into single pdf
outF_all <- sprintf("%s/plot_essCandidates_allFigs.pdf", outDir)
pdf(outF_all, width=10, height=5)
print(plot_list)
print(p1)
print(p2)
print(p4)
print(p5)
dev.off()
