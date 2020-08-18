# Script to test predictive capability of VBC score on guide effectiveness
library(openxlsx)
library(data.table)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(jtools)

####
## I/O
####

# VBC score flile
# https://www.nature.com/articles/s41592-020-0850-8#code-availability
vbc_file <- "input_data/vbc_score/hg38_all_sgRNAs.txt"

# Reduced TKOv3 library files
tko_file <- "input_data/TKOv3/newLibrary_tkov3_gRNAs_max_FINAL.xlsx"
tpam_file <- "input_data/TKOv3/gRNA_seqs30_pam.txt"

# Depmap LFC files
deplfc_file <- "input_data/depmap/Achilles_logfold_change_20Q2.csv"
dpam_file <- "input_data/depmap/gRNA_seqs30_pam.txt"

# Tiling LFC files
vlfc_files2 <- list.files(path = "input_data/val", pattern = "foldchange_mean", full.names = TRUE)
vpam_file <- "input_data/val/gRNA_seqs30_pam.txt"

# Core essentials
ess_file <- "input_data/essentials/core_essGenes.rda"
ness_file <- "input_data/essentials/nessGenes.rda"

# Output folder
output_folder <- "output_data/out_vbc"
if (!dir.exists(output_folder)) { dir.create(output_folder) }

####
## PREPARE DATA
####

# VBC (LOF prediction) data
vbc <- fread(vbc_file, data.table = FALSE)
colnames(vbc)[1] <- "GENE"
vbc2 <- vbc[,-15] # remove duplicated column
colnames(vbc2) <- gsub("-|=|\\ |>", "_", colnames(vbc2))
rm(vbc)
gc()

# Reduced TKOv3 library data
tko <- read.xlsx(tko_file)
tko2 <- tko[,c(2,3,6)]
colnames(tko2) <- c("GENE", "GUIDE", "LFC")

# Guide + PAM sequenec
tpam <- read.delim(tpam_file, h = FALSE, stringsAsFactors = FALSE)
tpam <- tpam[,c(1,3,4)]
colnames(tpam) <- c("GENE", "GUIDE", "sgRNA")
tpam$sgRNA <- substr(tpam$sgRNA, 5, nchar(tpam$sgRNA)-3)

# Merge with tko
tko2 <- left_join(tko2, tpam, by = c("GENE", "GUIDE"))

# Depmap LFC data (mean across cell lines)
deplfc <- fread(deplfc_file, data.table = FALSE)
deplfc <- data.frame(GUIDE = deplfc[,1], LFC = rowMeans(deplfc[,-1], na.rm = TRUE))

# Guide + PAM sequence
dpam <- read.delim(dpam_file, h = FALSE, stringsAsFactors = FALSE)
dpam <- dpam[,c(1,3,4)]
colnames(dpam) <- c("GENE", "GUIDE", "sgRNA")
dpam$sgRNA <- substr(dpam$sgRNA, 5, nchar(dpam$sgRNA)-3)

# Join with deplfc
deplfc <- left_join(deplfc, dpam, by = "GUIDE")
deplfc <- deplfc[,c(3,1,2,4)]

# Tiling Guide-level LFC data (only need WT)
vlfc_wt <- read.delim(grep("WT", vlfc_files2, value = TRUE), h = TRUE, stringsAsFactors = FALSE)
vlfc_wt <- vlfc_wt[,-which(colnames(vlfc_wt) %in% "Guide.ID")]
colnames(vlfc_wt)[3] <- "LFC"
vlfc_wt <- na.omit(vlfc_wt)

# Separate gene name from guide sequence
vlfc_wt <- separate(vlfc_wt, col = "GENE_CLONE", sep = "_", into = c("GENE", "GUIDE"))

# Guide + PAM sequence
vpam <- read.delim(vpam_file, h = FALSE, stringsAsFactors = FALSE)
vpam <- vpam[,c(1,3,4)]
colnames(vpam) <- c("GENE", "GUIDE", "sgRNA")
vpam$sgRNA <- substr(vpam$sgRNA, 5, nchar(vpam$sgRNA)-3)

# Join with vlfc_wt
vlfc_wt <- left_join(vlfc_wt, vpam, by = c("GENE", "GUIDE"))

# Core essential gene data
load(ess_file)
load(ness_file)

####
## ANALYSIS
####

# Define libraries
vlfc_wt$LIBRARY <- "Tiling"
tko2$LIBRARY <- "TKOv3_new"
deplfc$LIBRARY <- "Depmap"

# Combine data
dat_lfc <- bind_rows(vlfc_wt, tko2, deplfc)

# Indicate core essential and non-essential genes
dat_lfc$ESSENTIALITY <- "library_gene"
dat_lfc[which(dat_lfc$GENE %in% essGenes), "ESSENTIALITY"] <- "core_essential"
dat_lfc[which(dat_lfc$GENE %in% nessGenes), "ESSENTIALITY"] <- "non_essential"

# Original guide / gene numbers
n_org <- dat_lfc %>%
  group_by(LIBRARY) %>%
  summarise(genes = length(unique(GENE)), guides = n(), .groups = "keep")

# Combine library LFC with vbc score data
dat_vbc <- left_join(dat_lfc, vbc2, by = c("GENE", "sgRNA"))
dat_vbc2 <- dat_vbc[,c(1:6,10,12:14,20)]

# Remove data for guides without VBC score or LFC data
dat_vbc2 <- na.omit(dat_vbc2)
dat_vbc2 <- unique(dat_vbc2)

# Missing filtered guide / gene numbers
n_filt <- dat_vbc2 %>%
  group_by(LIBRARY) %>%
  summarise(genes = length(unique(GENE)), guides = n(), .groups = "keep")

# Remove data for genes with <= 2 guides
n_guides <- dat_vbc2 %>%
  group_by(LIBRARY, GENE) %>%
  summarise(N = n(), .groups = "keep") %>%
  as.data.frame

# Define genes with <= 2 guides
low_guides <- filter(n_guides, N <= 2)
dat_vbc2 <- anti_join(dat_vbc2, low_guides[,1:2], by = c("GENE", "LIBRARY"))

# Guide filtered guide / gene numbers
n_filt2 <- dat_vbc2 %>%
  group_by(LIBRARY) %>%
  summarise(genes = length(unique(GENE)), guides = n(), .groups = "keep")

# Generate randomized guide dropout data to compare
set.seed(42)
tko_rand <- subset(dat_vbc2, LIBRARY == "TKOv3_new")
tko_rand <- transform(tko_rand, LFC = sample(LFC))
tko_rand$LIBRARY <- "Randomized"

# Combine
dat_vbc3 <- bind_rows(dat_vbc2, tko_rand)

# Correlation between dropout and vbc, per gene
cor_gene <- dat_vbc3 %>%
  group_by(LIBRARY, GENE) %>%
  summarise(CORR = cor(LFC, VBC_score), N_GUIDE = n(), .groups = "keep") %>%
  arrange(CORR) %>%
  as.data.frame

####
## PLOT
####

# Set factor levels
dat_vbc3$ESSENTIALITY <- factor(
  dat_vbc3$ESSENTIALITY, levels = c("core_essential", "non_essential", "library_gene")
)

dat_vbc3$LIBRARY_F <- factor(
  dat_vbc3$LIBRARY, levels = c("Tiling", "TKOv3_new", "Depmap", "Randomized")
)

# Group colours
cols <- brewer.pal(3, "Dark2")
cols <- c(cols[1:2], "darkgrey")

# VBC score distribution
p1 <- ggplot(vbc2, aes(VBC_score)) +
        geom_density(fill = "grey") +
        geom_vline(aes(xintercept = median(VBC_score)), linetype = "dashed") +
        labs(title = "VBC score distribution", x = "VBC score", y = "Density",
             subtitle = sprintf("%s guides for %s genes", nrow(dat_vbc2), length(unique(dat_vbc2$GENE)))) +
        theme_minimal()

# Dropout vs vbc score
# Median LFC per group
lfc_median <- dat_vbc3 %>%
  group_by(LIBRARY_F) %>%
  summarise(median = median(LFC), N = n(), .groups = "keep") %>%
  as.data.frame

# Median LFC per group
vbc_median <- dat_vbc3 %>%
  group_by(LIBRARY_F) %>%
  summarise(median = median(VBC_score), N = n(), .groups = "keep") %>%
  as.data.frame

# Plot
p2 <- ggplot(dat_vbc3, aes(x = LFC, y = VBC_score, colour = ESSENTIALITY)) +
        facet_wrap(. ~ LIBRARY_F, scales = "free", ncol = 2, nrow = 2) +
        geom_point(alpha = 0.2) +
        geom_smooth(formula = y ~ x, method = "lm") +
        stat_ellipse() +
        stat_cor(method = "pearson", size = 3) +
        geom_hline(data = vbc_median, aes(yintercept = median), linetype = "dashed", colour = "black") +
        geom_vline(data = lfc_median, aes(xintercept = median), linetype = "dashed", colour = "black") +
        scale_colour_manual(values = cols) +
        labs(title = "Guide dropout vs. VBC score (LOF prediction) correlation",
             x = "Guide LFC", y = "VBC score",
             subtitle = sprintf("Tiling (%s guides, %s genes)\nTKOv3 (%s guides, %s genes)\nDepmap (%s guides, %s genes)",
                        nrow(dat_vbc3[which(dat_vbc3$LIBRARY == "Tiling"),]),
                        length(unique(dat_vbc3[which(dat_vbc3$LIBRARY == "Tiling"),"GENE"])),
                        nrow(dat_vbc3[which(dat_vbc3$LIBRARY == "TKOv3_new"),]),
                        length(unique(dat_vbc3[which(dat_vbc3$LIBRARY == "TKOv3_new"),"GENE"])),
                        nrow(dat_vbc3[which(dat_vbc3$LIBRARY == "Depmap"),]),
                        length(unique(dat_vbc3[which(dat_vbc3$LIBRARY == "Depmap"),"GENE"])))) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank())

# Print out
print(p2)

# Pick top 3 and worst 3 guides using vbc score and plot correlation with LFC
# demonstrate effectiveness of prediction score
val_vbc <- subset(dat_vbc3, LIBRARY == "Tiling")

# Rank guides per gene by LFC
val_vbc <- val_vbc %>%
  group_by(GENE) %>%
  arrange(GENE, LFC) %>%
  mutate(RANK = 1:n()) %>%
  as.data.frame

# Best 3 guides via vbc
val_best_vbc <- val_vbc %>%
  group_by(GENE) %>%
  arrange(GENE, desc(VBC_score)) %>%
  slice(1:3) %>%
  as.data.frame

# Worst 3 guides via vbc
val_worst_vbc <- val_vbc %>%
  group_by(GENE) %>%
  arrange(GENE, VBC_score) %>%
  slice(1:3) %>%
  as.data.frame

# Random 3 guides
set.seed(42)
val_rand_vbc <- val_vbc %>%
  group_by(GENE) %>%
  sample_n(3) %>%
  as.data.frame

# Combine to plot
val_best_vbc$LIBRARY <- "Tiling_best3_VBC"
val_worst_vbc$LIBRARY <- "Tiling_worst3_VBC"
val_rand_vbc$LIBRARY <- "Tiling_random3_VBC"
val_select <- bind_rows(val_best_vbc, val_worst_vbc, val_rand_vbc)

# Set factor levels
val_select$LIBRARY_F <- factor(
  val_select$LIBRARY,
  levels = c("Tiling_best3_VBC", "Tiling_worst3_VBC", "Tiling_random3_VBC")
)

# Median LFC per subset
lfc_median2 <- val_select %>%
  group_by(LIBRARY_F) %>%
  summarise(median = median(LFC), N = n(), .groups = "keep") %>%
  as.data.frame

# Median vbc per subset
vbc_median2 <- val_select %>%
  group_by(LIBRARY_F) %>%
  summarise(median = median(VBC_score), N = n(), .groups = "keep") %>%
  as.data.frame

# Plot scatter
p3 <- ggplot(val_select, aes(x = LFC, y = VBC_score, colour = ESSENTIALITY)) +
        facet_wrap(. ~ LIBRARY_F) +
        geom_point(alpha = 0.2) +
        geom_smooth(formula = y ~ x, method = "lm") +
        stat_ellipse() +
        stat_cor(method = "pearson", size = 3) +
        geom_hline(data = vbc_median2, aes(yintercept = median), linetype = "dashed", colour = "black") +
        geom_vline(data = lfc_median2, aes(xintercept = median), linetype = "dashed", colour = "black") +
        scale_colour_manual(values = cols) +
        labs(x = "Guide LFC", y = "VBC score") +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.title = element_blank())

# Plot eCDF
p4 <- ggplot(val_select, aes(RANK, colour = LIBRARY_F)) +
        stat_ecdf(geom = "smooth") +
        labs(x = "Guide rank", y = "Cumulative fraction") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title = element_blank())

# Plot density
p5 <- ggplot(val_select, aes(RANK, fill = LIBRARY_F)) +
        geom_density(alpha = 0.5) +
        labs(x = "Guide rank", y = "Density") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.title = element_blank())

# Merge plots and print out
p6 <- plot_grid(p4, p5, ncol = 2, rel_widths = c(1, 1))
p7 <- plot_grid(p3, p6, nrow = 2, rel_heights = c(1.3, 1))

# Print out
print(p7)

# Summarise model
# Fit linear regression model
fit <- lm(VBC_score ~ LFC, data = subset(dat_vbc3, LIBRARY == "Tiling"))
res <- fit$residuals
hist(res)
#plot(dat_vbc3$LFC, res)
#plot(dat_vbc3$VBC_score, res)

# Visualize model predictions
summ(fit)
effect_plot(fit, pred = LFC, interval = TRUE, plot.points = TRUE)

####
## REDUCED TKOv3 VS VBC
####

# Prepare data
tko3 <- tko2[,1:3]
colnames(tko3)[2] <- "SEQUENCE"

# Join tko and vbc data
vbc3 <- vbc2[,c(1:6)]
vbc3$SEQUENCE <- substr(vbc3$sgRNA, 0, nchar(vbc3$sgRNA)-3)
tko_vbc <- full_join(tko3, vbc3, by = c("GENE", "SEQUENCE"))
tko_vbc2 <- tko_vbc[,-c(19:23)]
tko_vbc2 <- tko_vbc2[which(tko_vbc2$GENE %in% tko3$GENE),]
tko_vbc2$guide_set <- ifelse(tko_vbc2$SEQUENCE %in% tko3$SEQUENCE, "TKOv3_new", "VBC")

# Order by gene and vbc
tko_vbc3 <- tko_vbc2 %>%
  group_by(GENE) %>%
  arrange(GENE, LFC, desc(VBC_score)) %>%
  mutate(mean_LFC = mean(LFC, na.rm = TRUE)) %>%
  as.data.frame

# Pick top VBC-scored guide for duplicated sequences
tko_vbc3 <- tko_vbc3[!duplicated(tko_vbc3$SEQUENCE),]

# Plot overall density distribution
p8 <- ggplot(tko_vbc3, aes(x = VBC_score, fill = guide_set)) +
        geom_density(alpha = 0.5, na.rm = TRUE) +
        theme_bw()

# Draw out
plot_out8 <- sprintf("%s/plot_density_TKOv3new_VBCscore.pdf", output_folder)
pdf(plot_out8, width = 8, height = 4, useDingbats = FALSE)
print(p8)
dev.off()

# VBC 99 quantile
vbc_quant <- as.numeric(quantile(tko_vbc3$VBC_score, 0.99, na.rm = TRUE))

# Qualify overall guide prediction
tko_vbc3[which(tko_vbc3$mean_LFC > -0.5 & tko_vbc3$VBC_score > vbc_quant & tko_vbc3$Auto_pick_top_sgRNAs <= 3), "guide_pred"] <- "weakLFC_highVBC"
tko_vbc3[which(tko_vbc3$mean_LFC < -0.5 & tko_vbc3$VBC_score > vbc_quant & tko_vbc3$Auto_pick_top_sgRNAs <= 3), "guide_pred"] <- "strongLFC_highVBC"
tko_vbc3[which(tko_vbc3$mean_LFC < -0.5 & tko_vbc3$VBC_score < vbc_quant & tko_vbc3$Auto_pick_top_sgRNAs > 3), "guide_pred"] <- "strongLFC_lowVBC"
tko_vbc3[which(tko_vbc3$mean_LFC > -0.5 & tko_vbc3$VBC_score < vbc_quant & tko_vbc3$Auto_pick_top_sgRNAs > 3), "guide_pred"] <- "weakLFC_lowVBC"

# Get number of genes in each group
to_plot <- tko_vbc3[,c("GENE", "guide_set", "guide_pred")]
to_plot <- na.omit(unique(to_plot))
to_plot <- as.data.frame(table(to_plot$guide_set, to_plot$guide_pred))

# Plot
p9 <- ggplot(to_plot, aes(x = Var1, y = Freq, fill = Var2)) +
        facet_grid(. ~ Var1, scales = "free") +
        geom_bar(stat = "identity", position = position_dodge()) +
        labs(x = NULL, y = "Number of genes") +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())

# Draw out
plot_out9 <- sprintf("%s/plot_bars_lowLFC_highVBC.pdf", output_folder)
pdf(plot_out9, width = 8, height = 4, useDingbats = FALSE)
print(p9)
dev.off()

# Filter for guides with low mean LFC and high VBC prediction
# For genes with multiple guides, sort by Auto_pick and select top one
tko_vbc4 <- tko_vbc3 %>%
  filter(guide_set == "VBC", guide_pred == "weakLFC_highVBC") %>%
  group_by(GENE) %>%
  arrange(Auto_pick_top_sgRNAs) %>%
  slice(1) %>%
  as.data.frame

# Fix tko column names
col_names <- colnames(tko)[-ncol(tko)]
rownames(tko) <- tko[,1]
tko <- tko[,-1]
colnames(tko) <- col_names

# Add vbc gudies to reduced TKOv3 guide file
vbc_df <- matrix(nrow = nrow(tko_vbc4), ncol = ncol(tko))
colnames(vbc_df) <- colnames(tko)
vbc_df <- as.data.frame(vbc_df)
vbc_df$GENE <- tko_vbc4$GENE
vbc_df$SEQUENCE <- tko_vbc4$SEQUENCE
vbc_df$Source <- "VBC predicted"

# Combine
rownames(tko) <- NULL
final <- bind_rows(tko, vbc_df)
final <- final[order(final$GENE),]
rownames(final) <- make.unique(final$GENE)

# Write out
table_out <- sprintf("%s/%s_VBC_CR.xlsx", output_folder, gsub(".xlsx", "", basename(tko_file)))
write.xlsx(table_final, file = table_out, row.names = FALSE)

# NOTE find incompatible guides (KB)
# The recipe I use for finding these is to pad them with CGACCG on the front
# and GTTTAG on the back, and then use a regular expression (or strict string
# matching will work) to search for CGTCTC or GAGACG in the full padded sequence.
# So you'll notice a bunch of the "failures" start with TCTC, which will match
# CGTCTC when dropped into the vector, and thus be cut incorrectly by the
# restriction enzymes during cloning.

blah <- final
blah$SEQUENCE_PAD <- paste("CGACCG", blah$SEQUENCE, "GTTTAG", sep = "")
grep("CGTCTC|GAGACG", blah$SEQUENCE_PAD)
