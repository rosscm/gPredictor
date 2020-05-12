# Script to plot CRISPR pipeline scores
library(plyr)
library(dplyr)
library(OneR)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggrepel)
library(openxlsx)
library(pheatmap)

dataDir <- "/Users/catherineross/GIN"
outDir <- sprintf("%s/bin/R/guideAnalysis", dataDir)

####################################################
# DATA INPUT
####################################################
# Output from prepOut.R
tko <- read.xlsx(sprintf("%s/output_data/out_TKOv3/table_TKOv3_guideScores.xlsx", outDir))
val <- read.xlsx(sprintf("%s/output_data/out_val/table_val_guideScores.xlsx", outDir))

## Merge data for plotting
## NOTE need to re-calculate gene fraction bins to make consistent across libraries
dat_merge <- rbind(tko, val)
dat_merge$gene_fraction_bin <- bin(dat_merge$gene_fraction, nbins=10, method="length")

####################################################
# PLOT
####################################################
# Various score boxplots
set=val
## CRISPRO scores
score_df <- set[,c("guide", "logFC", "essentiality", "doench_bin", "provean_bin", "disorder_bin", "guide_target", "SecStruct")]
score_df_melt <- melt(score_df, measure.vars=c("doench_bin", "provean_bin", "disorder_bin", "guide_target", "SecStruct"))
## Reorder value by factor levels
doench_lev <- levels(factor(score_df$doench_bin))
prov_lev <- levels(factor(score_df$provean_bin))
prov_lev <- prov_lev[c(5,3,2,1,7,6,4,8,9,10)]
dis_lev <- levels(factor(score_df$disorder_bin))
guide_lev <- levels(factor(score_df$guide_target))
struc_lev <- levels(factor(score_df$SecStruct))
score_df_melt$value <- factor(score_df_melt$value, levels=c(doench_lev, prov_lev, dis_lev, guide_lev, struc_lev))

## Guide positional info
score_df2 <- set[,c("guide", "logFC", "essentiality", "Exon", "gene_fraction_bin")]
score_df2 <- filter(score_df2, Exon <= 30) ## remove exon > 30
score_df_melt2 <- melt(score_df2, measure.vars=c("Exon", "gene_fraction_bin"))
## Reorder value by factor levels
exon_lev <- levels(factor(score_df2$Exon))
gene_lev <- levels(factor(score_df2$gene_fraction_bin))
score_df_melt2$value <- factor(score_df_melt2$value, levels=c(gene_lev, exon_lev))

## inDelphi, FORECasT, Azimuth scores
score_df3 <- set[,c("guide", "logFC", "essentiality", "frameshift_bin", "indel_bin", "efficiency_bin", "predicted_oof_bin")]
score_df_melt3 <- melt(score_df3, measure.vars=c("frameshift_bin", "indel_bin", "efficiency_bin", "predicted_oof_bin"))
## Reorder value by factor levels
frame_lev <- levels(factor(score_df3$frameshift_bin))
frame_lev <- frame_lev[c(7,1:6,8:10)]
indel_lev <- levels(factor(score_df3$indel_bin))
indel_lev <- indel_lev[c(3:5,1:2)]
eff_lev <- levels(factor(score_df3$efficiency_bin))
oof_lev <- levels(factor(score_df3$predicted_oof_bin))
oof_lev <- oof_lev[c(3,1:2,4:10)]
score_df_melt3$value <- factor(score_df_melt3$value, levels=c(frame_lev, indel_lev, eff_lev, oof_lev))

## Other
set$seqScore_bin <- bin(set$sequence_score, nbins=10, method="length")
score_df4 <- set[,c("guide", "logFC", "essentiality", "seqScore_bin")]
score_df_melt4 <- melt(score_df4, measure.vars="seqScore_bin")
## Reorder value by factor levels
seq_lev <- levels(factor(score_df4$seqScore_bin))
score_df_melt4$value <- factor(score_df_melt4$value, levels=seq_lev)
score_df4 <- set[,c("guide", "logFC", "essentiality", "seqScore_bin")]
score_df_melt4 <- melt(score_df4, measure.vars="seqScore_bin")
## Reorder value by factor levels
seq_lev <- levels(factor(score_df4$seqScore_bin))
score_df_melt4$value <- factor(score_df_melt4$value, levels=seq_lev)

## Plot
dat_plot <- score_df_melt3
p <- ggplot(dat_plot, aes(x=value, y=logFC, fill=essentiality)) +
        facet_grid(.~variable, scales="free_x", space="free_x") +
        geom_hline(yintercept=0, linetype="dashed") +
        geom_boxplot(outlier.alpha=0.3, outlier.size=1) +
        labs(x="Score bins", y="Mean guide-level LFC\n(across tiling screens)",
             title="Distribution of mean LFC across CRISPRO score bins") +
        theme_bw() +
        scale_fill_brewer(palette="Set2") +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
print(p)

plotF <- "~/Desktop/score.png"
ggsave(plotF, p, width=13, height=6)
