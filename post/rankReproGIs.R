# Script to prioritize guides targeting PDCD5/FASN interacting genes for query construction
library(openxlsx)
library(plyr)
library(dplyr)
library(OneR)
library(ggplot2)
library(gridExtra)

dataDir <- "/Users/catherineross/GIN"
inDir <- sprintf("%s/bin/R/guideAnalysis/output_data/out_TKOv3/plot_hypomorph", dataDir)
outDir <- inDir

#########################
## INPUT DATA
#########################
## Reproducible PDCD5/FASN interactions
datF <- sprintf("%s/PDCD5 and FASN reproducible interactions.xlsx", inDir)
dat <- read.xlsx(datF)

# PDCD5 interactions
dat1 <- dat[,c(1:8)]
colnames(dat1)[c(1,3:6)] <- c("query", "screen_1", "screen_2", "screen_3", "screen_4")
dat1 <- dat1[-which(is.na(dat1$query)),]

# FASN interactions
dat2 <- dat[,c(9:15)]
colnames(dat2)[c(1,3:5)] <- c("query", "screen_1", "screen_2", "screen_3")
dat2$screen_4 <- NA
dat2 <- dat2[,c(1:5,8,6,7)]

# Combine data
df <- rbind(dat1, dat2)

#########################
## GIN wildtype data
#########################
# Wt dropout (foldchange) data
#wtF <- list.files(pattern="GIN.*.xlsx", path=sprintf("%s/data/wildtypes", dataDir), full.names=TRUE)
#wt_name <- paste("wt", substr(basename(wtF), 4, 6), sep="_")
#wt <- lapply(wtF, function(x) read.xlsx(x, sheet=2))

#for (i in seq_along(wt)) {
#  tmp <- wt[[i]]
#  colnames(tmp)[4:ncol(tmp)] <- paste(wt_name[i], colnames(tmp)[4:ncol(tmp)], sep="_")
#  wt[[i]] <- tmp
#}

# Join data and keep relevant columns
#to_rm <- paste(c("Rich", "pyruvate", "PTC"), collapse="|")
#tps <- paste(c("T18", "T19", "T20"), collapse="|")
#wt_all <- join_all(wt, by=c("Guide.ID", "GENE_CLONE", "GENE"))
#wt_dat <- wt_all[,-grep(to_rm, colnames(wt_all))]
#wt_dat <- wt_dat[,grep(tps, colnames(wt_dat))]
#colnames(wt_dat) <- substr(colnames(wt_dat), 0, nchar(colnames(wt_dat))-12)

# Get mean wt log2FC
#wt_mean <- as.data.frame(rowMeans(wt_dat))
#wt_mean <- cbind(guide=wt_all$GENE_CLONE, gene=wt_all$GENE, wt_mean)
#colnames(wt_mean)[3] <- "logFC"
#wt_mean$guide <- gsub(".*\\_", "", wt_mean$guide) # get only guide sequence

load("wt_mean.rda")

# Get mean logFC per gene
wt_gene_mean <- aggregate(wt_mean$logFC, list(wt_mean$gene), mean, na.rm=TRUE)
wt_gene_mean <- na.omit(wt_gene_mean)
colnames(wt_gene_mean) <- c("gene", "mean_logFC")

#########################
## CRISPR guide score data
#########################
# Integrate with guide score data (crispro, azimuth, indelphi)
tkoF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_TKOv3/plot_hypomorph/table_TKOv3_guideScores.xlsx", dataDir)
valF <- sprintf("%s/bin/R/guideAnalysis/output_data/out_val/plot_hypomorph/table_tiling_guideScores.xlsx", dataDir)
tko <- read.xlsx(tkoF)
val <- read.xlsx(valF)

# All guide score data (TKOv3 + tiling)
dat_merge <- rbind(tko, val)
dat_merge$gene_fraction_bin <- bin(dat_merge$gene_fraction, nbins=10, method="length")
colnames(dat_merge)[c(2,66)] <- c("gene", "guide_efficiencyPred")

#########################
## Data aggregation
#########################
df_all <- join_all(list(df, wt_gene_mean, dat_merge), by="gene")
df_all$guide_library <- factor(df_all$guide_library, levels=c("TKOv3", "Tiling"))

#########################
## ANALYSIS
#########################
## Prioritize genes with consistent qGI score across screens
all_qGI_var <- c()
for (i in 1:length(unique(df_all$gene))) {
  df_gene <- filter(df_all, gene == unique(df_all$gene)[i])
  qGI_var <- var(c(df_gene[1,"screen_1"], df_gene[1,"screen_2"], df_gene[1,"screen_3"], df_gene[1,"screen_4"]), na.rm=TRUE)
  qGI_var <- round(qGI_var, 3)
  all_qGI_var[i] <- qGI_var
}

df_qGI <- data.frame(gene=unique(df_all$gene), qGI_var=all_qGI_var)
df_qGI <- df_qGI[order(df_qGI$qGI_var),]

## Interested in guides that target a protein domain
df_stat <- df_all[,c("query", "gene", "guide", "logFC", "sd", "dataset", "guide_library", "gene_fraction_bin", "Exon",
                     "guide_target", "Interpro_Description", "sequence_score", "guide_efficiencyPred")]
df_qGI_stat <- join(df_qGI, df_stat)
df_qGI_stat <- df_qGI_stat[,c(3,1,4,7,8,2,5,6,9:14)]
df_qGI_stat <- df_qGI_stat[with(df_qGI_stat, order(qGI_var, gene, logFC, -sequence_score, -guide_efficiencyPred)),]

## Filter table for guides targeting a domain, +logFC < -1.5, +positive sequence score
df_qGI_stat2 <- filter(df_qGI_stat, guide_target == "protein_domain")
df_qGI_stat3 <- filter(df_qGI_stat2, logFC <= -1.5)
df_qGI_stat4 <- filter(df_qGI_stat3, sequence_score >= 0)

## Write out
df_list <- list(df_qGI_stat, df_qGI_stat2, df_qGI_stat3, df_qGI_stat4)
outF <- sprintf("%s/table_TKOv3_tiling_reproducibleGIs_prioritized.xlsx", outDir)
write.xlsx(df_list, outF)

#########################
## PLOT
#########################
## Mean LFC vs. gene position (exon / gene fraction)
plot_list <- list()
df <- df_qGI_stat
df$guide_library <- factor(df$guide_library, levels=c("TKOv3", "Tiling"))

for (k in unique(df$gene)) {
  print(k)
  dat_plot <- filter(df, gene == k)

  if (nrow(dat_plot) == 0) next

  # Plot title + x-axis
  query <- unique(dat_plot$query)
  title <- sprintf("%s (%s reproducible GI)", k, query)
  stat <- "Exon"

  # Get gene-level mean LFC in each dataset
  means = ddply(dat_plot, ~dataset, summarise, mean=mean(logFC))
  dat_plot$geneLFC <- NA
  dat_plot[which(dat_plot$dataset == "TKOv3"),]$geneLFC <- means[which(means$dataset == "TKOv3"),]$mean
  dat_plot[which(dat_plot$dataset == "Tiling"),]$geneLFC <- means[which(means$dataset == "Tiling"),]$mean

  # Manual colour scale
  cols <- c("#7570b3", "#1b9e77")
  names(cols) <- levels(df$guide_library)
  colScale <- scale_colour_manual(name="guide_library", values=cols)

  plot_list[[k]] <-
    ggplot(dat_plot, aes(x=factor(Exon), y=logFC, colour=guide_library, shape=Interpro_Description)) +
        facet_wrap(.~dataset) +
        geom_point(alpha=0.5, aes(size=sequence_score)) +
        geom_errorbar(aes(ymin=logFC-sd, ymax=logFC+sd), width=0.2, alpha=0.5) +
        colScale +
        scale_size_continuous(range=c(1,4)) +
        geom_hline(aes(yintercept=geneLFC), linetype="dashed", colour="blue") +
        labs(x=stat, title=title,
             y="Mean guide-level LFC\n(across min WT or tiling screens)",
             subtitle=sprintf("Guide LFC distribution vs %s (y=gene LFC)", stat)) +
        theme_bw() +
        theme(text=element_text(family="sans", size=15))
}

outF <- sprintf("%s/point_TKOv3_tiling_reproducibleGIs_prioritized_guideLFC_vs_%s_perGene.pdf", outDir, stat)
plot_grob <- marrangeGrob(plot_list, nrow=5, ncol=3)
ggsave(outF, plot_grob, width=40, height=25)
