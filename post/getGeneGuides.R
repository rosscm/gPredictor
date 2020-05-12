# Get list of guides for set of genes (for Elvin, Dick lab)
library(openxlsx)
library(data.table)
library(tidyverse)

######
# INPUT
######

# List of query genes
genes_file <- "input_data/dick_lab/30 genes_Elvin_DickLab.xlsx"
genes <- read.xlsx(genes_file)
genes <- c(genes[,1], colnames(genes), "OR2B11")

## NOTE OR2W5 (or OR2W5P) not found in crispro annotations
## Using OR2B11 instead - enhancer gene for OR2W5

# Sequence score file (Hart lab)
seq_score_file <- "input_data/gRNA_seqScores_Hart_2017.txt" # gRNA sequence score table (BAGEL 2017 paper)
seq_score <- read.delim(seq_score_file, h = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

# Crispro annotation file
anno_file <- "software/crispro-master/annotations/CRISPRO.GRCh37.ensembl90.20180901.csv"
anno <- fread(anno_file, h = TRUE, data.table = FALSE)

######
# GUIDE LIST PER GENE
######
######

# Subset crispro annotations for query genes
genes2 <- paste(paste0("^", genes, "$"), collapse = "|")
gene_anno <- anno[grep(genes2, anno$gene_name),]

# Remove entries without guide sequence
gene_anno2 <- gene_anno[which(is.na(gene_anno$guide) == FALSE),]

# Filter out MINOR/ALTERNATIVE isoforms (APPRIS column)
gene_anno3 <- gene_anno2[grep("PRINCIPAL", gene_anno2$APPRIS),]

# Remove duplicated guides
gene_anno4 <- gene_anno3[!duplicated(gene_anno3$guide),]

# Only keep relevant columns
gene_anno4 <- gene_anno4[,-c(6:9,15,27:29,33,34,36)]

######
# CALCULATE GUIDE SEQUENCE SCORE
######

# Modify sequence score table for easier indexing
rownames(seq_score) <- seq_score$X
seq_score$X <- NULL
colnames(seq_score) <- seq_along(seq_score)

# Calculate sequence score for all guides
guides <- gene_anno4$guide
guide_sequence <- strsplit(guides, "") # list of guides split by nucleotide

guide_scores <- c()
for (i in seq_along(guide_sequence)) {
  print(i)
  guide <- guide_sequence[[i]]
  score <- c()
  for (position in seq_along(guide)) {
    nucleotide <- guide[position]
    score <- sum(c(score, seq_score[nucleotide, position]))
    guide_scores[i] <- score
  }
}

guide_scores_df <- data.frame(
  guide = guides,
  gene_name = gene_anno4$gene_name,
  sequence_score = guide_scores
)

# Merge sequence scores with annotations
dat <- left_join(gene_anno4, guide_scores_df)

######
# TKOv3
######

# File of TKOv3 guides
tko_guide_file <- "/Users/catherineross/GIN/data/pipeline_input/experimentalInfo/TKOv3_Library_20170518.txt"
tko_guide <- read.delim(tko_guide_file, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

# Join with dat
dat$TKOv3 <- "No"
dat[which(dat$guide %in% tko_guide$SEQUENCE), "TKOv3"] <- "Yes"

######
# GIN WT LFC
######

# GIN Wildtype LFC data
wt_file <- list.files(pattern="GIN.*.xlsx", path="/Users/catherineross/GIN/data/wildtypes", full.names=TRUE)
wt_name <- paste("wt", substr(basename(wt_file), 4, 6), sep = "_")
wt <- lapply(wt_file, function(x) read.xlsx(x, sheet = 2))

# Rename columns to indiciate analysis number for merging
for (i in seq_along(wt)) {
  tmp <- wt[[i]]
  colnames(tmp)[4:ncol(tmp)] <- paste(wt_name[i], colnames(tmp)[4:ncol(tmp)], sep="_")
  wt[[i]] <- tmp
}

# Join wt data together and keep relevant columns
#to_rm <- paste(c("Rich", "pyruvate", "PTC"), collapse="|")
to_rm <- paste(c("pyruvate", "PTC"), collapse = "|") # keep rich columns
tps <- paste(c("T18", "T19", "T20"), collapse = "|")
wt_all <- reduce(wt, left_join, by = c("Guide.ID", "GENE_CLONE", "GENE"))
wt_dat <- wt_all[,-grep(to_rm, colnames(wt_all), ignore.case = TRUE)]
wt_dat <- wt_dat[,grep(tps, colnames(wt_dat), ignore.case = TRUE)]

# Get mean wt log2FC, variance, and sd
### Minimal
wt_min <- wt_dat[grep("Minimal", colnames(wt_dat), ignore.case = TRUE)]
wt_mean_min <- rowMeans(wt_min, na.rm = TRUE)
wt_var_min <- apply(wt_min, 1, var, na.rm = TRUE)
wt_sd_min <- apply(wt_min, 1, sd, na.rm = TRUE)

### Rich
wt_rich <- wt_dat[grep("Rich", colnames(wt_dat), ignore.case = TRUE)]
wt_mean_rich <- rowMeans(wt_rich, na.rm = TRUE)
wt_var_rich <- apply(wt_rich, 1, var, na.rm = TRUE)
wt_sd_rich <- apply(wt_rich, 1, sd, na.rm = TRUE)

### Combine both wt sets
wt_mean <- data.frame(
  guide = wt_all$GENE_CLONE, gene_name = wt_all$GENE,
  logFC_min = wt_mean_min, var_min = wt_var_min, sd_min = wt_sd_min,
  logFC_rich = wt_mean_rich, var_rich = wt_var_rich, sd_rich = wt_sd_rich
)
wt_mean$guide <- gsub(".*\\_", "", wt_mean$guide) # get only guide sequence

# Join with dat
dat <- left_join(dat, wt_mean, by = c("gene_name", "guide"))

######
# PREPARE FINAL DATA
######

# Sort by sequence score per gene to get top-ranked gRNA per gene
dat2 <- dat %>%
    group_by(gene_name) %>%
    arrange(gene_name, desc(sequence_score)) %>%
    as.data.frame()

# Generate list of dfs per gene
dat_list <- list()
for (i in unique(dat2$gene_name)) {
  dat_i <- filter(dat2, gene_name == i)
  dat_list[[i]] <- dat_i
}

# Write out excel workbook
fname <- gsub("\\..*", "", basename(genes_file))
fname <- gsub(" ", "_", fname)
table_out <- sprintf("table_%s_gRNAs_CR.xlsx", fname)
write.xlsx(dat_list, file = table_out, row.names = FALSE)
