# Script to plot guide scores and test for hypomorph-generating guides
# Load all packages
packages <- c("dplyr", "OneR", "reshape2", "ggplot2", "gridExtra"
              "ggpubr", "ggrepel", "openxlsx")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Set data directory
output_folder <- "~/projects/gPredictor/output_data"

##############
# DATA INPUT
##############

# Output from prepOut.R
tko_file <- sprintf("%s/out_TKOv3/table_TKOv3_guideScores_whole.xlsx", output_folder)
val_file <- sprintf("%s/out_val/table_val_guideScores_whole.xlsx", output_folder)
tko <- read.xlsx(tko_file)
val <- read.xlsx(val_file)

# Merge data for plotting
## NOTE need to re-calculate gene fraction bins to make consistent across libraries

## NOTE from 06-28-2021 ##
## tko and val have different logFC columns (tko calculated for min/rich)
## val only calculated for one condition
## rbind call no longer works (fix if need to use code again ...)

dat_merge <- rbind(tko, val)
dat_merge$gene_fraction_bin <- bin(dat_merge$gene_fraction, nbins=10, method="length")

##############
# DATA ANALYSIS
##############

# Compare same guide in TKOv3 and tiling library
tko_fc <- tko[,c("guide", "gene_name", "logFC")]
colnames(tko_fc)[3] <- paste(colnames(tko_fc)[3], "TKO", sep="_")
val_fc <- val[,c("guide", "gene_name", "logFC")]
colnames(val_fc)[3] <- paste(colnames(val_fc)[3], "tiling", sep="_")

# Combine logFC data for guides in both TKOv3 and tiling and get difference
fc_comb <- join(tko_fc, val_fc)
fc_comb <- na.omit(fc_comb)
fc_comb$fc_diff <- abs(fc_comb$logFC_TKO - fc_comb$logFC_tiling)
fc_comb <- fc_comb[order(fc_comb$fc_diff, decreasing=TRUE),]
fc_comb$int <- 1:nrow(fc_comb)

# Label guides >1.5 foldchange difference
fc_comb$label <- NA
toLabel <- which(fc_comb$fc_diff >= 1.5)
fc_comb[toLabel,]$label <- paste(fc_comb[toLabel,]$gene_name, fc_comb[toLabel,]$guide, sep="_")

## Plot difference in LFC measurements across libraries
out_file <- sprintf("%s/point_TKOv3_tiling_LFC_consistency.pdf", output_folder)
pdf(out_file, width=13, height=5)
p <- ggplot(fc_comb, aes(x=int, y=fc_diff, label=label)) +
        geom_point(size=1.5, alpha=0.5) +
        labs(x="Guide (in both libraries)", y="Diff LFC",
             title=sprintf("LFC consistency across TKOv3 and tiling libraries (N=%s guides)", nrow(fc_comb))) +
        geom_hline(aes(yintercept=mean(fc_diff)), linetype="dashed", colour="blue") +
        geom_text_repel(parse=TRUE, na.rm=TRUE, segment.size=0.2, segment.color="grey50",
                        nudge_x=500-subset(fc_comb, label != "NA")$int,
                        size=3, direction="y", hjust=0) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank())
print(p)
dev.off()

##############
# GENERATE PLOTS
##############

## 1) mean log2FC vs. exon / gene fraction
plotPoint <- function(data, name, ess=TRUE, stats) {
  plot_list_1 <- list()

  if (ess==TRUE) {
    dat_filt <- filter(data, essentiality == "essential")
    out_file_1 <- sprintf("%s/point_%s_ess_guideLFC_vs_%s_perGene.pdf", output_folder, name, stats)
  }
  if (ess==FALSE) {
    dat_filt <- data
    out_file_1 <- sprintf("%s/point_%s_all_guideLFC_vs_%s_perGene.pdf", output_folder, name, stats)
  }

  for (k in unique(dat_filt$gene_name)) {
    print(k)
    dat_plot <- filter(dat_filt, gene_name == k)
    # annotate whether or not gene is core essential
    if (k %in% core$V1) {
      title <- sprintf("%s (core essential)", k)
    } else {
      title <- sprintf("%s (not core essential)", k)
    }

    plot_list_1[[k]] <-
      ggplot(dat_plot, aes(x=factor(get(stats)), y=logFC, colour=guide_library, shape=guide_target)) +
        geom_point(size=2, alpha=0.5) +
        geom_errorbar(aes(ymin=logFC-sd, ymax=logFC+sd), width=0.2, alpha=0.5) +
        scale_colour_manual(values=c("black", "red")) +
        geom_hline(aes(yintercept=meanLFC), linetype="dashed", colour="blue") +
        labs(x=stats, y="Mean guide-level LFC\nacross queries", title=title, subtitle=sprintf("Guide LFC distribution vs %s", stats)) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
  }
  cat("Arranging plots together...")
  plot_grob_1 <- marrangeGrob(plot_list_1, nrow=5, ncol=5)
  ggsave(out_file_1, plot_grob_1, width=40, height=35)
  cat(" done.")
}

plotPoint(tko, name="wt", ess=TRUE, stats="gene_fraction_bin")
plotPoint(tko, name="wt", ess=FALSE, stats="gene_fraction_bin")

## 2) same as #1 but with TKO + VAL data combined (for CB)
plot_list_2 <- list()
for (k in unique(dat_merge$gene_name)) {
  print(k)
  #k = "URB1"  #"DDB1" "URB1" "RPS19BP1" (spot checking genes)
  dat_plot <- filter(dat_merge, gene_name == k)
  stat <- "gene_fraction"

  if (unique(dat_plot$essentiality) == "essential") {
    title <- sprintf("%s (essential,", k)
    if (unique(dat_plot$core) == "Core") {
      title <- sprintf("%s core)", title)
    }
    if (unique(dat_plot$core) == "Not_core") {
      title <- sprintf("%s not core)", title)
    }
  } else {
    title <- sprintf("%s (non-essential)", k)
  }

  # Get gene-level mean LFC in each dataset
  means = ddply(dat_plot, ~dataset, summarise, mean=mean(logFC))
  dat_plot$geneLFC <- NA
  dat_plot[which(dat_plot$dataset == "TKOv3"),]$geneLFC <- means[which(means$dataset == "TKOv3"),]$mean
  dat_plot[which(dat_plot$dataset == "Tiling"),]$geneLFC <- means[which(means$dataset == "Tiling"),]$mean

  plot_list_2[[k]] <-
    ggplot(dat_plot, aes(x=gene_fraction_bin, y=logFC, colour=guide_library, shape=guide_target)) +
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
out_file_2 <- sprintf("%s/point_TKOv3_tiling_all_guideLFC_vs_%s_perGene.pdf", output_folder, stat)
plot_grob_2 <- marrangeGrob(plot_list_2, nrow=5, ncol=3)
ggsave(out_file_2, plot_grob_2, width=40, height=35)

## 3) rank plots of all guides in library + strip chart
# plot
out_file_3 <- sprintf("%s/strip_guideLFC_perExon_perGene.pdf", output_folder)
pdf(out_file_3, width=20, height=15, useDingbats=FALSE)
#for (k in unique(val_ess$gene_name)) {
#  print(k)
#  # Removing data for exons with <100 guides targeting it
#  dat_plot <- filter(val_ess, Exon <= 30 & gene_name == k)

#  p3 <- ggplot(val_ess, aes(x=Rank, y=logFC)) +
#          facet_wrap(.~Exon, ncol=6) +
#          geom_point(alpha=0.3) +
#          geom_rug(data=dat_plot, sides="b", aes(colour=guide_library)) +
#          scale_colour_manual(values=c("black", "green")) +
#          geom_hline(yintercept=0, linetype="dashed", colour="blue") +
#          labs(x="Guide rank", y="Mean LFC", title=k,
#               subtitle="Guide LFC distribution across library, per exon") +
#          theme_bw() +
#          theme(text=element_text(family="sans", size=15))
#  print(p3)
#}
#dev.off()

## 4) distribution of gene fraction vs logFC
# point
out_file_4 <- sprintf("%s/point_guideLFC_gene_fraction", output_folder)
pdf(sprintf("%s.pdf", out_file_4), width=15, height=7.5)
p4 <- ggplot(val_ess, aes(x=gene_fraction, y=logFC, colour=guide_target)) +
        geom_point(alpha=0.25) +
        geom_smooth(se=TRUE) +
        labs(x="Gene fraction", y="Mean LFC", title="Gene fraction vs. meanLFC") +
        scale_colour_manual(values=c("black", "lightblue")) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15))
print(p4)
dev.off()

# facet by chromosome
out_file_4b <- sprintf("%s_perChrom.pdf", out_file_4)
pdf(out_file_4b, width=20, height=15)
p4b <- p4 + facet_wrap(.~chrom, ncol=5)
print(p4b)
dev.off()

# binned boxplot
out_file_5 <- sprintf("%s/boxplot_guideLFC_vs_gene_fraction", output_folder)
pdf(sprintf("%s.pdf", out_file_5), width=5.5, height=7)
p5 <- ggplot(val_ess, aes(x=gene_fraction_bin, y=logFC)) +
        geom_boxplot(fill="lightgrey") +
        geom_hline(yintercept=0, linetype="dashed", colour="blue") +
        labs(x="Gene fraction", y="Mean LFC", title="Gene fraction vs. meanLFC") +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
print(p5)
dev.off()

# facet by chromosome
out_file_5b <- sprintf("%s_perChrom.pdf", out_file_5)
pdf(out_file_5b, width=15, height=15)
p5b <- p5 + facet_wrap(.~chrom, ncol=5)
print(p5b)
dev.off()

## 6) fraction of guide LFC < 1.5 per exon
df <- data.frame()
for (j in 1:max(val_ess$Exon)) {
  df[j,1] <- j
  dat_exon <- filter(val_ess, Exon == j)
  if (!length(dat_exon)) {
    df[j,2] <- 0
  } else {
    df[j,2] <- round((length(which(dat_exon$logFC < -1.5))/nrow(dat_exon))*100, 2)
  }
}
colnames(df) <- c("Exon", "Fraction")

out_file_7 <- sprintf("%s/bar_guideFraction_vs_Exon.pdf", output_folder)
pdf(out_file_7, width=30, height=5)
p7 <- ggplot(df, aes(x=factor(Exon), y=Fraction)) +
        geom_bar(stat="identity") +
        labs(x="Exon", y="Guide fraction",
             title=sprintf("Fraction of guides with strong dropout (LFC < -1.5) across %s exons", max(val_ess$Exon))) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15))
print(p7)
dev.off()

## 7) additional score plots (eg disorder, provean, oof, etc.)
# Mean crispro scores across gene fraction bins
means <- aggregate(val[,c("doench_score", "provean_score", "disorder_score")], list(val$gene_fraction_bin, val$essentiality), mean)
means_melt <- melt(means)

out_file_8 <- sprintf("%s/bar_guideScores_vs_gene_fraction.pdf", output_folder)
pdf(out_file_8, width=17, height=4.5)
p8 <- ggplot(means_melt, aes(x=Group.1, y=value)) +
        facet_wrap(.~variable, scales="free", nrow=1) +
        geom_bar(stat="identity") +
        labs(x="Gene fraction", y="Score mean",
             title="Mean of guide-level scores across gene position") +
        theme_bw() +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
print(p8)
dev.off()

## 8) mean LFC across crispro score bins
# boxplot
# CRISPRO scores
score_df <- val[,c("guide", "logFC", "essentiality", "doench_bin", "provean_bin", "disorder_bin", "guide_target", "SecStruct")]
score_df_melt <- melt(score_df, measure.vars=c("doench_bin", "provean_bin", "disorder_bin", "guide_target", "SecStruct"))

out_file_9 <- sprintf("%s/boxplot_guideLFC_vs_crispro_scores.png", output_folder)
pdf(out_file_9, width=20, height=5)
p9 <- ggplot(score_df_melt, aes(x=factor(value), y=logFC, fill=essentiality)) +
        facet_grid(.~variable, scales="free_x", space="free_x") +
        geom_hline(yintercept=0, linetype="dashed") +
        geom_boxplot(outlier.alpha=0.3, outlier.size=1) +
        labs(x="Score bins", y="Mean guide-level LFC\n(across tiling screens)",
             title="Distribution of mean LFC across CRISPRO score bins") +
        theme_bw() +
        scale_fill_brewer(palette="Set2") +
        theme(text=element_text(family="sans", size=15),
              axis.text.x=element_text(angle=90, hjust=1))
print(p9)
dev.off()

# point
score_df2 <- val_ess[,c("guide", "logFC", "guide_target", "doench_score", "oof_score", "disorder_score", "provean_score")]
score_df_melt2 <- melt(score_df2, measure.vars=c("doench_score", "oof_score", "disorder_score", "provean_score"))

out_file_10 <- sprintf("%s/point_guideLFC_vs_crispro_scores.pdf", output_folder)
pdf(out_file_10, width=20, height=5)
p10 <- ggplot(score_df_melt2, aes(x=value, y=logFC, colour=guide_target)) +
        facet_wrap(.~variable, scales="free", nrow=1) +
        geom_point(alpha=0.25) +
        geom_smooth(se=TRUE) +
        labs(x="Score", y="Mean LFC", title="Guide-level mean LFC across CRISPRO scores") +
        scale_colour_manual(values=c("black", "lightblue")) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15))
print(p10)
dev.off()

## 9) mean LFC per protien secondary structure
#dat_plot <- filter(val, SecStruct != "None")
dat_plot <- val

out_file_11 <- sprintf("%s/boxplot_guideLFC_vs_SecStruct.pdf", output_folder)
pdf(out_file_11, width=15, height=7)
p11 <- ggplot(dat_plot, aes(x=SecStruct, y=logFC)) +
        facet_wrap(.~essentiality) +
        geom_boxplot() +
        labs(x="Secondary structure prediction", y="Mean LFC", title="Guide-level mean LFC across protein\nsecondary structure predictions") +
        theme_bw() +
        theme(text=element_text(family="sans", size=15)) +
        stat_compare_means(label="p.signif", method="t.test", ref.group="None")
print(p11)
dev.off()
