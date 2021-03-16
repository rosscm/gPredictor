# Script to filter gRNA library by CRISPR scores
#' @param lib (char) name of library being tested (ie TKOv3, val/tiling)

library(plyr)
library(dplyr)
library(openxlsx)
library(pheatmap)
library(ggplot2)

dataDir <- "/Users/catherineross/GIN"

# Function to filter for most effective guides
predGuides <- function(lib, essential=TRUE, tko_guides=FALSE,
                      seqScore_filt=0, azmScore_filt=0.6,
                      doenchScore_filt=60) {

  # Generate output directory
  dt <- format(Sys.Date(), "20%y%m%d")
  outDir  <- sprintf("%s/res/%s_out_%s_predGuides", dataDir, dt, lib)
  if (!file.exists(outDir)) dir.create(outDir)

  # Read in score table
  datF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_%s/table_%s_guideScores.xlsx", dataDir, lib, lib)
  dat <- read.xlsx(datF)

  # Run on essential genes only or whole set (TRUE->essential only)
  if (essential == TRUE) dat <- filter(dat, essentiality == "essential")
  # Keep or filter TKOv3 guides from tiling set (TRUE->keep)
  if (tko_guides == FALSE) dat <- filter(dat, guide_library != "TKOv3")

  # Generate matrix to test correlation between logFC and CRISPR scores
  # ie which scores to use to filter gRNA library
  mat <- dat[,c("logFC", "Exon", "gene_fraction", "doench_score", "provean_score",
                "disorder_score", "Frameshift.frequency",
                "Expected.indel.length", "azimuth_guideEfficiency_pred",
                "sequence_score", "predicted_oof")]
  name <- paste(dat$gene_name, dat$guide, sep="_")
  rownames(mat) <- name

  # Measure correlation
  score_cor <- cor(mat, use="pairwise.complete.obs")
  hc <- hclust(as.dist(1-score_cor))
  plot_mat <- score_cor[hc$order, hc$order]
  plot_mat[upper.tri(plot_mat)] = NA

  # Plot heatmap
  plotF1 <- sprintf("%s/plot_corr_CRISPRscores_logFC_%s", outDir, lib)
  if (essential == TRUE) plotF1 <- sprintf("%s_essential", plotF1)

  png(sprintf("%s.png", plotF1), width=8, height=8, units="in", res=600)
  print(pheatmap(plot_mat, display_numbers=TRUE, cluster_col=FALSE, cluster_row=FALSE,
      main="Correlation between log2FC and CRISPR gRNA scores"))
  dev.off()

  # Scores most 'predictive' of guide dropout: doench, efficiency, sequence score
  ## Select best guides per gene based on scores
  guide_list <- c()
  guide_df <- c()
  for (i in unique(dat$gene_name)) {
    gene <- dat[which(dat$gene_name == i),]
    cat(sprintf("Filtering guides for gene: %s\n", i))
    cat(sprintf("  - total guides: %s\n", nrow(gene)))

    if (nrow(gene) <= 1) {
      # Get guide ID and write out to object
      guides <- paste(gene$gene_name, gene$guide, sep="_")
      guide_list <- c(guide_list, guides)
      guide_df <- rbind(guide_df, gene)
    } else {
      # Filter by sequence score, keep positives
      gene_filt <- filter(gene, sequence_score > seqScore_filt)
      cat(sprintf("  - filtering for sequence score > %s: %s\n", seqScore_filt, nrow(gene_filt)))
      # Filter by guide efficiency prediction >= 60%
      gene_filt <- filter(gene_filt, azimuth_guideEfficiency_pred >= azmScore_filt)
      cat(sprintf("  - filtering for efficiency prediction >= %s: %s\n", azmScore_filt, nrow(gene_filt)))
      # Filter by doench score >= 60
      gene_filt <- filter(gene_filt, doench_score >= doenchScore_filt)
      cat(sprintf("  - filtering for doench score >= %s: %s\n\n", doenchScore_filt, nrow(gene_filt)))
      # Filter by protein domain target
      #gene_filt <- filter(gene_filt, guide_target == "protein_domain")
      #cat(sprintf("  - filtering for protein domain targets: %s\n\n", nrow(gene_filt)))

      # Get guide IDs for remaining guides and write out to object
      guides <- paste(gene_filt$gene_name, gene_filt$guide, sep="_")
      guide_list <- c(guide_list, guides)
      guide_df <- rbind(guide_df, gene_filt)
    }
  }

  #filt <- filter(dat,
  #  gene_fraction <= 0.5 &
  #  doench_score >= doenchScore_filt &
  #  provean_score < -5 &
  #  disorder_score < 0.5 &
  #  guide_target == "protein_domain" &
  #  azimuth_guideEfficiency_pred >= azmScore_filt &
  #  Frameshift.frequency >= 50 &
  #  SecStruct != "None" &
  #  sequence_score > seqScore_filt
  #)

  # Compare filtered guide distribution with whole
  df_filt <- guide_df[,c("gene_name", "guide", "logFC", "essentiality"),]
  df_filt$set <- "Scoring_subset"
  df_all <- dat[,c("gene_name", "guide", "logFC", "essentiality"),]
  df_all$set <- "Total_guides"
  df <- rbind(df_filt, df_all)

  # Stats
  cat(sprintf("No. total guides: %s; genes: %s\n", nrow(df_all), length(unique(df_all$gene_name))))
  cat(sprintf("No. predicted subset guides: %s; genes: %s\n", nrow(df_filt), length(unique(df_filt$gene_name))))

  #blah <- filt[which(filt$logFC >= 0),] # check which guides still have weak dropout after filtering
  pval <- ks.test(x=df_filt$logFC, y=df_all$logFC, alternative="greater")
  print(pval)

  # Mean of filtered and un-filtered guide sets
  means <- ddply(df, "set", summarise, grp.mean=mean(logFC))

  # Plot distributions
  p2 <- ggplot(df, aes(x=logFC, fill=set)) +
          geom_density(alpha=0.5) +
          labs(y="Density", x=sprintf("Mean guide-level LFC (across %s screens)", lib),
               title="Distibution of whole vs. filtered guide set",
               subtitle="Filtered by sequence score, predicted efficiency, doench score") +
          geom_vline(data=means, aes(xintercept=grp.mean, color=set), linetype="dashed", size=1) +
          theme_bw() +
          scale_fill_brewer(palette="Set1") +
          theme(text=element_text(family="sans", size=14),
                legend.position="bottom")

  plotF2 <- sprintf("%s/plot_density_CRISPRscores_logFC_%s", outDir, lib)
  if (essential == TRUE) plotF2 <- sprintf("%s_essential", plotF2)
  ggsave(sprintf("%s.png", plotF2), p2, width=6.5, height=4)

  # Return vector of guides that meet set criteria above
  return(guide_list)
}

# Run function
#val_guides <- predGuides(lib="val")
predGuides(lib="val")
