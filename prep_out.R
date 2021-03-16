# Script to write out table with all CRISPR tool guide scores for specified library

packages <- c("plyr", "OneR", "openxlsx", "data.table", "argparse")
for (p in packages) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

######
# USER INPUTS
######

guide_library = "val"

# User files
essentials_file <- "~/data/essentials/DepMap_essential_19Q2_60_percent.txt"
core_file       <- "~/data/essentials/ess_traver.txt" # core essentials
tko_guide_file  <- "~/data/TKOv3/TKOv3_Library_20170518.txt" # TKOv3 guides
seq_score_file <- "input_data/TKOv3/gRNA_seqScores_Hart_2017.txt" # gRNA sequence score table (BAGEL 2017 paper)
vbc_file       <- "input_data/vbc_score/hg38_all_sgRNAs.txt"

# gPredictor files
input_file     <- sprintf("input_data/%s/gRNA_guideSeq.txt", guide_library) # original guide library file
crispro_file   <- sprintf("output_data/out_%s/crispro/CRISPRO.GRCh37.ensembl90.20180901_%s.txt", guide_library, guide_library)
indel_file     <- sprintf("output_data/out_%s/indelphi/indelphi_K562_guideSeq79_stats.txt", guide_library)
azimuth_file   <- sprintf("output_data/out_%s/azimuth/azimuth_guideSeq30_predictions.txt", guide_library)
forecast_file  <- sprintf("output_data/out_%s/forecast/forecast_all_predictedindelsummary.txt", guide_library)

# Output directory
output_folder <- sprintf("output_data/out_%s", guide_library)

######
# READ IN FILES
######

# Read in files
essentials <- read.delim(essentials_file, h = FALSE, as.is = TRUE)
core <- read.delim(core_file, h = FALSE, as.is = TRUE)
tko_guide <- read.delim(tko_guide_file, h = TRUE, as.is = TRUE)
seq_score <- read.delim(seq_score_file, h = TRUE, as.is = TRUE)
vbc <- fread(vbc_file, h = TRUE, data.table = FALSE)
input <- read.delim(input_file, h = TRUE, as.is = TRUE)
crispro <- read.delim(crispro_file, h = TRUE, as.is = TRUE)
indel <- read.delim(indel_file, h = TRUE, as.is = TRUE)
azimuth <- read.delim(azimuth_file, h = TRUE, as.is = TRUE)
forecast <- fread(forecast_file, fill = TRUE, data.table = FALSE)

######
# MODIFY TABLES FOR MERGING
######

# Guide sequence scores
# Modify sequence score table for easier indexing
rownames(seq_score) <- seq_score$X
seq_score$X <- NULL
colnames(seq_score) <- seq_along(seq_score)

# CRISPRO
# NOTE need to deal with multiple gene transcripts
# First, filter for principal isoform as indicated by APRIS column
crispro2 <- crispro[order(crispro$APPRIS, decreasing=TRUE),]
# Remove irrelevant columns
crispro2 <- crispro2[,-c(33,34,36)]
crispro3 <- crispro2[!duplicated(crispro2$guide),] # remove duplicated guides

# inDelphi
indel <- indel[,-3]
colnames(indel)[1:2] <- c("gene_name", "guide")

# Azimuth
azimuth <- azimuth[,-3]
colnames(azimuth) <- c("gene_name", "guide", "guideEfficiency_pred")

# FORECAST
# Get guide info rows (first column preceded by "@@@")
# @@@id guide_seq predicted_in_frame
forecast_guide <- forecast[grep("@@@", forecast$V1),]
forecast_guide <- unique(forecast_guide)
colnames(forecast_guide) <- c("guide", "predicted_in_frame", "predicted_oof")
forecast_guide$predicted_in_frame <- as.numeric(forecast_guide$predicted_in_frame)
# Get original guide IDs from inDelphi output (retained original order)
forecast_guide$guide <- indel$guide
# Calculate out-of-frame prediction (100-(predicted-in-frame) per guide)
forecast_guide$predicted_oof <- 100-(forecast_guide$predicted_in_frame)

# VBC (guide efficacy prediction)
# Grab relevant columns
vbc2 <- vbc[,c(1,2,6,7)]
colnames(vbc2) <- c("gene_name", "guide", "vbc_score", "distance_tss_0_stop_1")
# Remove PAM sequence
vbc2$guide <- substr(vbc2$guide, 0, nchar(vbc2$guide)-3)

######
# GET FOLDCHANGE DATA (TKOv3 or TILING)
######

# Dropout (foldchange) data
## TKOv3
if (guide_library == "TKOv3") {

  wt_file <- list.files(pattern="GIN.*.xlsx", path="~/data/wildtypes", full.names=TRUE)
  wt_name <- paste("wt", substr(basename(wt_file), 4, 6), sep="_")
  wt <- lapply(wt_file, function(x) read.xlsx(x, sheet=2))

  # Rename columns to indiciate analysis number for merging
  for (i in seq_along(wt)) {
    tmp <- wt[[i]]
    colnames(tmp)[4:ncol(tmp)] <- paste(wt_name[i], colnames(tmp)[4:ncol(tmp)], sep="_")
    wt[[i]] <- tmp
  }

  # Join wt data together and keep relevant columns
  #to_rm <- paste(c("Rich", "pyruvate", "PTC"), collapse="|")
  to_rm <- paste(c("pyruvate", "PTC"), collapse="|") # keep rich columns
  tps <- paste(c("T18", "T19", "T20"), collapse="|")
  wt_all <- join_all(wt, by=c("Guide.ID", "GENE_CLONE", "GENE"))
  wt_dat <- wt_all[,-grep(to_rm, colnames(wt_all), ignore.case=TRUE)]
  wt_dat <- wt_dat[,grep(tps, colnames(wt_dat), ignore.case=TRUE)]

  # Get mean wt log2FC, variance, and sd
  ### Minimal
  wt_min <- wt_dat[grep("Minimal", colnames(wt_dat), ignore.case=TRUE)]
  wt_mean_min <- rowMeans(wt_min, na.rm=TRUE)
  wt_var_min <- apply(wt_min, 1, var, na.rm=TRUE)
  wt_sd_min <- apply(wt_min, 1, sd, na.rm=TRUE)

  ### Rich
  wt_rich <- wt_dat[grep("Rich", colnames(wt_dat), ignore.case=TRUE)]
  wt_mean_rich <- rowMeans(wt_rich, na.rm=TRUE)
  wt_var_rich <- apply(wt_rich, 1, var, na.rm=TRUE)
  wt_sd_rich <- apply(wt_rich, 1, sd, na.rm=TRUE)

  ### Combine both wt sets
  wt_mean <- data.frame(
    guide=wt_all$GENE_CLONE, gene_name=wt_all$GENE,
    logFC_min=wt_mean_min, var_min=wt_var_min, sd_min=wt_sd_min,
    logFC_rich=wt_mean_rich, var_rich=wt_var_rich, sd_rich=wt_sd_rich
  )
  wt_mean$guide <- gsub(".*\\_", "", wt_mean$guide) # get only guide sequence
}

## Tiling
if (guide_library == "val") {

  wt_file <- list.files(pattern="foldchange_mean.txt", path="~/data/validation_tiling", full.names=TRUE)
  wt_name <- substr(basename(wt_file), 0, nchar(basename(wt_file))-29)
  wt <- lapply(wt_file, function(x) read.delim(x, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE))

  ## Merge query foldchange data together
  ## (expected essential genes behave similarly across different query backgrounds)
  wt_all <- join_all(wt, by=c("Guide.ID", "GENE_CLONE", "GENE"))

  ## Remove entries without guide info (05-11-2020)
  wt_all <- wt_all[-which(wt_all$GENE_CLONE == ""),]

  ## Prepare for getting mean values
  wt_dat <- wt_all[,-c(1:3)]
  names(wt_dat) <- wt_name

  # Get mean wt log2FC, variance, and sd
  wt_mean <- rowMeans(wt_dat, na.rm=TRUE)
  wt_var <- apply(wt_dat, 1, var, na.rm=TRUE)
  wt_sd <- apply(wt_dat, 1, sd, na.rm=TRUE)
  wt_mean <- data.frame(
    guide=wt_all$GENE_CLONE, gene_name=wt_all$GENE,
    logFC=wt_mean, var=wt_var, sd=wt_sd
  )
  wt_mean$guide <- gsub(".*\\_", "", wt_mean$guide) # get only guide sequence
}

######
# CALCULATE GUIDE SEQUENCE SCORE
######

# All library guides
guides <- input$SEQUENCE
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

guide_scores_df <- data.frame(guide=input$SEQUENCE, gene_name=input$GENE, sequence_score=guide_scores)

######
# CONSTRUCT FINAL SCORE TABLE
######

# Rename orginal input table columns
colnames(input) <- c("chrom", "start", "stop", "strand", "guide", "gene_name")

# Add additional library info
input <- cbind(
  input,
  dataset=guide_library,
  essentiality="non_essential",
  guide_library=guide_library,
  core="not_core"
)
input$essentiality <- as.character(input$essentiality)
input$guide_library <- as.character(input$guide_library)
input$core <- as.character(input$core)
input[which(input$gene_name %in% essentials$V1), "essentiality"] <- "essential"
input[which(input$guide %in% tko_guide$SEQUENCE), "guide_library"] <- "TKOv3"
input[which(input$gene_name %in% core$V1), "core"] <- "core"

# List of all guide feature tables
all_features <- list(input, crispro3, indel, azimuth, forecast_guide, guide_scores_df, vbc2)
if (exists("wt_mean")) {
  all_features[[length(all_features)+1]] <- wt_mean
}

# Join everything together
dat <- join_all(all_features)
dat <- unique(dat)

# Clean up domain target column
dat[is.na(dat$Interpro_Description), "Interpro_Description"] <- "None"
dat$guide_target <- "no_protein_domain"
dat[which(dat$Interpro_Description != "None"), "guide_target"] <- "protein_domain"

# Clean up secondary structure column
dat[is.na(dat$SecStruct), "SecStruct"] <- "None"
dat$guide_structure <- "no_protein_structure"
dat[which(dat$SecStruct != "None"), "guide_structure"] <- "protein_structure"

# Rank guides by strength of dropout (strongest=1)
if (guide_library == "TKOv3" | guide_library == "val") {
  if (guide_library == "TKOv3") { dat <- dat[order(dat$logFC_min),] }
  if (guide_library == "val")   { dat <- dat[order(dat$logFC),] }
  # Rank and revert to original order
  dat$Rank <- 1:nrow(dat)
  dat <- dat[order(as.numeric(rownames(dat))),]
}

# Get mean wt LFC of essentials and non-essentials
#means = ddply(dat, ~essentiality, summarise, mean=mean(logFC_min, na.rm=TRUE))
#dat$meanLFC_ess <- means[1,2]
#dat$meanLFC_ness <- means[2,2]

# Bin continuous data/score columns for easier plotting
dat$gene_fraction_bin <- bin(dat$gene_fraction, nbins=10, method="length", na.omit=FALSE)
dat$doench_bin <- bin(dat$doench_score, nbins=10, method="length", na.omit=FALSE)
dat$oof_bin <- bin(dat$oof_score, nbins=10, method="length", na.omit=FALSE)
dat$disorder_bin <- bin(dat$disorder_score, nbins=10, method="length", na.omit=FALSE)
dat$provean_bin <- bin(dat$provean_score, nbins=10, method="length", na.omit=FALSE)
dat$frameshift_bin <- (dat$Frameshift.frequency, nbins=10, method="length", na.omit=FALSE)
dat$indel_bin <- bin(dat$Expected.indel.length, nbins=5, method="length", na.omit=FALSE)
dat$efficiency_bin <- bin(dat$guideEfficiency_pred, nbins=10, method="length", na.omit=FALSE)
dat$predicted_oof_bin <- bin(dat$predicted_oof, nbins=10, method="length", na.omit=FALSE)
dat$sequence_score_bin <- bin(dat$sequence_score, nbins=10, method="length", na.omit=FALSE)
dat$vbc_score_bin <- bin(dat$vbc_score, nbins=10, method="length", na.omit=FALSE)

######
# WRITE OUT TABLE
######

fname <- sprintf("output_data/out_%s/table_%s_guideScores_whole", guide_library, guide_library)
write.xlsx(dat, file = paste0(fname, ".xlsx"))
write.table(dat, file = paste0(fname, ".txt"), quote = FALSE, sep = "\t")
