# Script to run enrichment on essential gene scores
library(plyr)
library(dplyr)
library(fgsea)
library(data.table)
library(reshape2)
library(ggplot2)
library(scales)

dataDir <- "/Users/catherineross/GIN"
inDir <- sprintf("%s/res/20191216_out_tiling_score_essGIs_predGuides", dataDir)
annoDir <- sprintf("%s/anno", dataDir)
outDir <- sprintf("%s/gsea", inDir)
if (!file.exists(outDir)) dir.create(outDir)

####################################################
# DATA INPUT
####################################################
## Pathway annotation set
annoF <- sprintf("%s/Human_GOBP_AllPathways_no_GO_iea_October_01_2018_symbol.gmt", annoDir)
## Input score data (essential gene LFC from tiling array)
scoreF <- sprintf("%s/table_tiling_scores_LFC_essGenes.txt", inDir)

####################################################
# DATA ANALYSIS
####################################################
# Run GSEA
runGSEA <- function(scoreF, annoF,
                    minGene=10, maxGene=500,
                    sigPath=FALSE, sigPathFDR=0.2) {

  # Read in pathway file
  path <- gmtPathways(annoF)
  # Filter by size
  path2 <- path[which(lapply(path, length) >= minGene)]
  path2 <- path2[which(lapply(path2, length)  <= maxGene)]

  # Prepare ranked data
  score <- read.delim(scoreF, h=TRUE, stringsAsFactors=FALSE)
  cat(sprintf("\nRead score data for %i interactions across %i query screens\n", nrow(score), ncol(score)))

  # Run GSEA per screen
  for (i in 1:ncol(score)) {

    cat(sprintf("\n[%i] Running GSEA on query screen %s...", i, names(score)[i]))

    # Create score vector and rank
    dat <- setNames(score[,i], rownames(score))
    datRnk <- sort(dat, decreasing=TRUE)
    # Remove NAs to avoid issue with ties in ranked list (2018-10-05)
    datRnk <- datRnk[!is.na(datRnk)]

    # Run GSEA (preranked)
    set.seed(42)
    gseaRun <- fgsea(stats=datRnk, pathways=path, nperm=10000, nproc=1)
    nOrgPath <- nrow(gseaRun)

    # Filtering for significant pathways (2018-11-19)
    if (sigPath == TRUE) {
      gseaRun <- as.data.table(filter(gseaRun, padj <= sigPathFDR))
      nSigPath <- nrow(gseaRun)
      cat(sprintf("\n**Filtered to %i enriched pathways from %i total (FDR <= %g)",
        nSigPath, nOrgPath, sigPathFDR))
    }

    nSigPath <- nrow(gseaRun)
    f_name <- substr(basename(scoreF), 0, nchar(basename(scoreF))-4)
    outF <- sprintf("%s/gsea_%s_%s", outDir, f_name, names(score)[i])

    # Split out positive / negative enrichments
    enrich_pos <- gseaRun[NES > 0][head(order(NES, decreasing=TRUE), n=30), pathway]
    enrich_neg <- gseaRun[NES < 0][head(order(NES, decreasing=FALSE), n=30), pathway]

    # Make GSEA table plots for negatively and positively enriched pathways
    if (nrow(gseaRun) == 0) {
      cat("\nNo significant pathways remaining!\n\n")
      next
    }

    cat(sprintf("\n  => plotting top positively enriched pathways to %s_pos_enrich.pdf.\n", basename(outF)))
    pdf(sprintf("%s_pos_enrich.pdf", outF), width=18, height=13)
    plotGseaTable(pathways=path[enrich_pos],
                  stats=datRnk,
                  fgseaRes=gseaRun,
                  gseaParam=0.5)
    invisible(dev.off())

    cat(sprintf("  => plotting top negatively enriched pathways to %s_neg_enrich.pdf.\n", basename(outF)))
    pdf(sprintf("%s_neg_enrich.pdf", outF), width=18, height=13)
    plotGseaTable(pathways=path[enrich_neg],
                  stats=datRnk,
                  fgseaRes=gseaRun,
                  gseaParam=0.5)
    invisible(dev.off())

    # Order results by +/- NES and write to output directory
    gseaRun_pos <- gseaRun[NES > 0][order(NES, decreasing=TRUE)]
    cat(sprintf("  => writing out positive enrichment results to %s_pos.txt.\n", basename(outF)))
    fwrite(gseaRun_pos, file=sprintf("%s_pos.txt", outF), sep="\t", sep2=c("", " ", ""))

    gseaRun_neg <- gseaRun[NES < 0][order(NES, decreasing=FALSE), ]
    cat(sprintf("  => writing out negative enrichment results to %s_neg.txt.\n", basename(outF)))
    fwrite(gseaRun_neg, file=sprintf("%s_neg.txt", outF), sep="\t", sep2=c("", " ", ""))

    # Plot gene-level heat maps
    plotHeatmap <- function(resSign) {

      if (resSign == "pos") { gseaRes <- gseaRun_pos }
      if (resSign == "neg") { gseaRes <- gseaRun_neg }

      if (nrow(gseaRes) == 0) {
        cat(sprintf("\nNo %s pathways to plot\n\n", resSign))
        return(NULL)
      }

      # Only plot top 30 pathways
      gsea_dat <- gseaRes$leadingEdge
      names(gsea_dat) <- gseaRes$pathway
      gsea_dat <- head(gsea_dat, n=30)

      # Melt to dataframe for ggplot
      gsea_dat_plot <- reshape2::melt(gsea_dat)
      colnames(gsea_dat_plot) <- c("Gene", "Pathway")

      # Transform pi scores vector to df and merge with gsea result df
      datRnk_df <- as.data.frame(datRnk)
      datRnk_df$Gene <- rownames(datRnk_df)
      rownames(datRnk_df) <- c()
      colnames(datRnk_df)[1] <- "LFC_score"

      # Prevent merge from re-arranging columns using plyr::join (2018-11-08)
      final <- join(gsea_dat_plot, datRnk_df, by="Gene")

      # Prevent ggplot from re-arranging those levels
      final$Pathway <- factor(final$Pathway, levels=gseaRes$pathway)
      final$Gene <- factor(final$Gene, levels=unique(final$Gene))

      # Write out heatmaps
      plotF <- sprintf("%s_%s_heatmap", outF, resSign)
      cat(sprintf("  => plotting LFC score heatmap to %s.pdf\n", basename(plotF)))

      ggplot(final, aes(Gene, Pathway)) +
          geom_tile(aes(fill=LFC_score), colour="white") +
          scale_fill_gradient2(low=muted("#386cb0"), high="#ef3b2c") +
          theme(panel.background=element_rect(fill="white"),
                panel.grid=element_blank(),
                plot.title=element_text(size=15, face="bold", hjust=0.5),
                axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8, face="bold"),
                axis.text.y=element_text(size=8, hjust=1, vjust=1, face="bold"),
                axis.title=element_text(size=12),
                axis.ticks=element_blank(),
                legend.title=element_text(face="bold", size=10)) +
          ggtitle(basename(plotF)) +
          labs(x="Leading edge gene", y="Pathway")

      ggsave(sprintf("%s.pdf", plotF), width=50, height=10, limitsize=FALSE)

     }
     plotHeatmap(resSign="pos") # heat maps for positively enriched pathways
     plotHeatmap(resSign="neg") # heat maps for negatively enriched pathways
  }
}

# Run function
runGSEA(scoreF, annoF, minGene=10, maxGene=500, sigPath=FALSE)
