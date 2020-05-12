#!/usr/bin/env Rscript

# Script to prepare input for run_gPredictor.R pipeline
# Input file needs following headers: CHROMOSOME START STOP STRAND SEQUENCE GENE
input_file <- "input_data/depmap/Achilles_guide_map_19Q4.csv"

if (grepl("csv", input_file)) {
  guides <- read.csv(input_file, stringsAsFactors=FALSE)
} else if (grepl("txt", input_file)) {
  guides <- read.delim(input_file, stringsAsFactors=FALSE)
} else if (grepl("xlsx", input_file)) {
  library(openxlsx)
  guides <- read.xlsx(input_file)
} else {
  cat("Unrecognized file format (supports txt, csv, xlsx formats).")
}

###########
# Reformat so data has the following headers: CHROMOSOME START STOP STRAND SEQUENCE GENE
## NOTE these steps need to be done on a case-by-case basis
# Extract relevant info from columns
guide_seq <- guides$sgrna
guide_gene <- unlist(lapply(strsplit(guides$gene, " "), "[", 1))
guide_align <- strsplit(guides$genome_alignment, "_")
guide_chr <- unlist(lapply(guide_align, "[", 1))
guide_strand <- unlist(lapply(guide_align, "[", 3))
guide_coord <- as.numeric(unlist(lapply(guide_align, "[", 2)))

# Put together final dataframe and write out
df <- data.frame(CHROMOSOME=guide_chr,
                 START=guide_coord,
                 STOP=guide_coord,
                 STRAND=guide_strand,
                 SEQUENCE=guide_seq,
                 GENE=guide_gene)

## NOTE using tiling data to help decipher proper depmap guide coordinates
# Negative strand coordinates
df[which(df$STRAND == "-"), "START"] <- df[which(df$STRAND == "-"), "START"] - 3
df[which(df$STRAND == "-"), "STOP"] <- df[which(df$STRAND == "-"), "START"] + unique(nchar(guide_seq))

# Positive strand coordinates
df[which(df$STRAND == "+"), "STOP"] <- df[which(df$STRAND == "+"), "STOP"] + 3
df[which(df$STRAND == "+"), "START"] <- df[which(df$STRAND == "+"), "STOP"] - unique(nchar(guide_seq))

###########

output_file <- sprintf("%s/gRNA_guideSeq.txt", dirname(input_file))
write.table(df, output_file, col=TRUE, row=FALSE, quote=FALSE, sep="\t")
