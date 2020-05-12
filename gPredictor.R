#!/usr/bin/env Rscript

# Script to analyze CRISPR guides using suite of softwares as outlined by He et al. 2019
# https://www.nature.com/articles/s41467-019-12489-8.pdf

### ====================================================
# Software descriptions
## Azimuth (https://github.com/MicrosoftResearch/Azimuth)
###    On-target guide prediction
## inDelphi (https://github.com/maxwshen/inDelphi-model)
###    Predicts DNA repair outcomes at double strand breaks
## CRISPRO (https://gitlab.com/bauerlab/crispro)
###   Maps functional scores of tiling sgRNAs to genomes, transcripts, protein
###   coordinates, and structures, providing general views of structure-function
###   relationships at discrete protein regions
## FORECasT or SelfTarget (https://github.com/felicityallen/SelfTarget)
###    Predicts CRISPR/Cas9-generated mutations
## ProTiler (https://github.com/MDhewei/ProTiler-1.0.0)
###   Maps protein regions that are associated with strong sgRNA dropout effect
###   in screens termed CRISPR knockout hyper-sensitive (CKHS) regions

###### ====================================================
######
# SET UP R ENVIRONMENT
######

pkgs <- c("reticulate", "breakfast", "stringr", "plyr", "dplyr", "data.table", "argparse")
for (p in pkgs) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Loads scoring functions
source(file.path("src", "functions.R"))

######
# PARSES USER ARGUMENTS
######

# Makes argparse object and arguments
parser <- ArgumentParser(description="Runs CRISPR guide scoring pipeline on specified set of gRNAs.")
parser$add_argument("-l", "--guide_library", type="character",
                    help="Name of guide library being scored (e.g., TKOv3) for output file naming")
parser$add_argument("-b", "--genome_build", type="character",
                    help="Genome build of guide library (hg19 or hg38)")
parser$add_argument("-i", "--input_file", type="character",
                    help="Path to file containing guide library information")
parser$add_argument("-o", "--output_folder", type="character", default="output_data",
                    help="Path to output folder [default %(default)s]")
args <- parser$parse_args()

######
# PARAMETER SETTING
######

guide_library <- args$guide_library
genome_build <- args$genome_build
input_file <- args$input_file
output_folder <- args$output_folder

# Checks if files exist
if (!file.exists(input_file)) {
  stop("ERROR: Input guide library file does not exist!")
}

# Install genome annotation correponding to specific genome build
if (genome_build == "hg19") genome_pkg <- "BSgenome.Hsapiens.UCSC.hg19"
if (genome_build == "hg38") genome_pkg <- "BSgenome.Hsapiens.NCBI.GRCh38"
suppressPackageStartupMessages(library(genome_pkg, character.only = TRUE))

# Create library-specific output folder if it doesn't already exist
output_folder <- "output_data"
if (!file.exists(output_folder)) dir.create(output_folder)
output_folder <- sprintf("%s/out_%s", output_folder, guide_library)
if (!file.exists(output_folder)) dir.create(output_folder)

######
# INSTALL SOFTWARE & DEPENDENCIES
######

## Call install.sh to install system dependencies and CRISPR analysis packages
## NOTE: using python2 for Azimuth compatibility
software_folder <- "software"
if (!file.exists(software_folder)) dir.create(software_folder)
#system("sh install.sh")

###### ====================================================
######
# EXTRACT GUIDE SEQUENCE CONTEXT
######

cat("### PREPARING INPUT DATA ###\n\n")

# Input file needs following headers: CHROMOSOME START STOP STRAND SEQUENCE GENE
cat(sprintf("* Reading input guide library file: %s\n", input_file))
guides <- read.delim(input_file, h=TRUE, as.is=TRUE, stringsAsFactors=FALSE)

cat(sprintf("** Guide library: %s\n", guide_library))
cat(sprintf("** Guide library build: %s\n", genome_build))
cat(sprintf("** Total guides in set: %s\n", nrow(guides)))
cat(sprintf("** Total targeted genes: %s\n\n", length(unique(guides$GENE))))

######
# FETCH 30BP GUIDE SEQUENCE CONTEXT
######

cat("* Extracting 30bp guide sequence context...\n")

# Define input parameters
input_list <- list(guides, guides)
build_list <- rep(genome_build, 2)
seq_list <- rep(30, 2)
strand_list <- c("+", "-")

seq30 <- mapply(getSeqContext,
  input=input_list,
  genome_build=build_list,
  seq_length=seq_list,
  strand=strand_list
)

seq30_df <- rbind(as.data.frame(seq30[,1]), as.data.frame(seq30[,2]))
seq30_df <- seq30_df[match(guides$SEQUENCE, seq30_df$guide),] # match original data order

# Stop code from continuing if sequence length is incorrect
seq30_len <- unique(unlist(lapply(as.character(seq30_df$x), nchar)))
if (seq30_len != 30) {
  stop("Guide sequence context extraction did not run properly!")
}

# Remove guides without proper 'NGG' PAM sequence
seq30_pam <- substr(as.character(seq30_df$pam), 2, nchar(as.character(seq30_df$pam)))
seq30_remove <- which(seq30_pam != "GG")

if (length(seq30_remove)) {
  warning(sprintf("Removing %s guides without proper 'NGG' PAM sequence\n", length(seq30_remove)))
  seq30_df <- seq30_df[-seq30_remove,]
}

seq30_file <- sprintf("%s/gRNA_seqs30_pam.txt", dirname(input_file))
write.table(seq30_df, seq30_file, col=FALSE, row=FALSE, quote=FALSE, sep="\t")

######
# FETCH 79BP GUIDE SEQUENCE CONTEXT
######

cat("* Extracting 79bp guide sequence context...\n")

# Define input parameters
input_list <- list(guides, guides)
build_list <- rep(genome_build, 2)
seq_list <- rep(79, 2)
strand_list <- c("+", "-")

seq79 <- mapply(getSeqContext,
  input=input_list,
  genome_build=build_list,
  seq_length=seq_list,
  strand=strand_list
)

seq79_df <- rbind(as.data.frame(seq79[,1]), as.data.frame(seq79[,2]))
seq79_df <- seq79_df[match(guides$SEQUENCE, seq79_df$guide),] # match original data order

# Stop code from continuing if sequence length / PAM sequence are incorrect
seq79_len <- unique(unlist(lapply(as.character(seq79_df$x), nchar)))
if (seq79_len != 79) {
  stop("Guide sequence context extraction did not run properly!")
}

# Remove guides without proper 'NGG' PAM sequence
seq79_pam <- substr(as.character(seq79_df$pam), 2, nchar(as.character(seq79_df$pam)))
seq79_remove <- which(seq79_pam != "GG")

if (length(seq79_remove)) {
  warning(sprintf("Removing %s guides without proper 'NGG' PAM sequence\n", length(seq79_remove)))
  seq79_df <- seq79_df[-seq79_remove,]
}

seq79_file <- sprintf("%s/gRNA_seqs79_pam.txt", dirname(input_file))
write.table(seq79_df, seq79_file, col=FALSE, row=FALSE, quote=FALSE, sep="\t")

###### ====================================================
######
# BEGIN ANALYSIS
######

cat("### INITIATING CRISPR GUIDE SCORING PIPELINE ###\n\n")

# Guide sequences and target genes (after processing)
guides <- guides[which(guides$SEQUENCE %in% seq30_df[,3]),]
seqs <- guides$SEQUENCE
genes <- guides$GENE

# Guide sequences + context (30bp and 79bp)
seqs30  <- as.character(seq30_df[,4])
seqs79  <- as.character(seq79_df[,4])

cat(sprintf("** Total guides after processing: %s\n", length(seqs)))
cat(sprintf("** Total targeted genes: %s\n", length(unique(genes))))

######
# AZIMUTH
######

message("\n\nRunning Azimuth for on-target guide efficiency prediction\n\n")

use_virtualenv("virtualenvs/azimuth", required = TRUE)
repl_python()
import azimuth.model_comparison
import numpy as np

seq = np.array(r.seqs30) # model trained with 30bp guides
seq_pred = azimuth.model_comparison.predict(seq, None, None)

# Return to R environment
exit
azimuth <- data.frame(GENE=genes, GUIDE=seqs, SEQUENCE=py$seq, PREDICTION=py$seq_pred)

# Write out results
azimuth_folder <- sprintf("%s/azimuth", output_folder)
azimuth_output <- sprintf("%s/azimuth_guideSeq30_predictions.txt", azimuth_folder)
if (!file.exists(azimuth_folder)) dir.create(azimuth_folder)
write.table(azimuth, azimuth_output, col=TRUE, row=FALSE, quote=FALSE, sep="\t")


######
# INDELPHI
######

message("\n\nRunning inDelphi to measure in-frame probabilities\n\n")

# Set up environment to import the inDelphi.py script
# Supported cell types are ['mESC', 'U2OS', 'HEK293', 'HCT116', 'K562']
celltype <- "K562"

use_virtualenv("virtualenvs/indelphi", required = TRUE)
repl_python()
import sys
sys.path.append("software/inDelphi-model")
import inDelphi
inDelphi.init_model(celltype=r.celltype)

# Iterate over gRNA sequences and store outputs in 'indel_pred' and 'indel_stats'
seq = r.seqs79
cutsite = 40 # 3 bases "before" PAM (0-based index)
indel_pred = []
indel_stats = []
for i in range(len(seq)):
  print(seq[i])
  pred_df, stats = inDelphi.predict(seq[i], cutsite)
  indel_pred.append(pred_df)
  indel_stats.append(stats)

# Return to R environment and prepare output
exit
indel_stats <- py$indel_stats
indel_stats_df <- lapply(indel_stats, function(x) as.data.frame(t(unlist(x))))
indel_stats_df <- as.data.frame(rbindlist(indel_stats_df))

# Re-order columns
indel_stats_df <- indel_stats_df[,c(3,1,2,4:ncol(indel_stats_df)),]
colnames(indel_stats_df)[1] <- "SEQUENCE"
indel_stats_df <- cbind(GENE=genes, GUIDE=seqs, indel_stats_df)

# Write out results
indel_folder <- sprintf("%s/indelphi", output_folder)
indel_output <- sprintf("%s/indelphi_%s_guideSeq79_stats.txt", indel_folder, celltype)
if (!file.exists(indel_folder)) dir.create(indel_folder)
write.table(indel_stats_df, indel_output, col=TRUE, row=FALSE, quote=FALSE, sep="\t")


######
# CRISPRO
######

message("\n\nRunning CRISPRO for extracting guide annotation information\n\n")

# NOTE cross-referencing input gRNA library with CRISPRO annotation file
# rather than running CRISPRO tool (does not run to completion)
# WARNING very large file (14GB)
crispro_file <- sprintf("%s/crispro-master/annotations/CRISPRO.GRCh37.ensembl90.20180901.csv", software_folder)
crispro <- fread(crispro_file, h=TRUE, data.table=FALSE)

# Merge with input data
guides_merge <- guides[,c("CHROMOSOME", "SEQUENCE", "GENE")]
colnames(guides_merge) <- c("chrom", "guide", "gene_name")
guides_crispro <- join(guides_merge, crispro)

# Write out results
crispro_folder <- sprintf("%s/crispro", output_folder)
crispro_output <- sprintf("%s/%s_%s.txt", crispro_folder, substr(basename(crispro_file), 0, nchar(basename(crispro_file))-4), guide_library)
if (!file.exists(crispro_folder)) dir.create(crispro_folder)
write.table(guides_crispro, crispro_output, col=TRUE, row=FALSE, quote=FALSE, sep="\t")


######
# FORECAST
######

message("\n\nRunning FORECAST to predict CRISPR/Cas9-generated mutations\n\n")

# Prepare batch file to run in Docker container
# NOTE PAM index is 0-based
guide_id  <- paste("Guide", seq(seqs79), genes, sep="_")
forecast_batch <- data.frame("ID"=guide_id, "Target"=seqs79, "PAM Index"=43, check.names=FALSE)
forecast_batch_file <- sprintf("%s/forecast_batch.txt", dirname(input_file))
write.table(forecast_batch, forecast_batch_file, col=TRUE, row=FALSE, quote=FALSE, sep="\t")

# Run in container
docker_call <- sprintf("sh forecast_docker.sh %s", forecast_batch_file)
system(docker_call)

# Read in output files
## NOTE finish once latest iteration completes (2020-02-04)
forecast_files <- list.files(pattern=".*batch.*indelsummary", path=sprintf("output_data/out_%s/forecast", guide_library), full.names=TRUE)

# Concatenate files and write out
forecast_files <- paste(forecast_files, collapse=" ")
forecast_output <- sprintf("output_data/out_%s/forecast/forecast_all_predictedindelsummary.txt", guide_library)
cat_call <- sprintf("cat %s >> %s", forecast_files, forecast_output)
system(cat_call)

### ====================================================
## END OF CODE
### ====================================================
