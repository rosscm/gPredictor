# Script to predict hypomorphic phenotypes from CRISPR guide data
library(data.table)
library(dplyr)
library(tidyr)
library(rio)
library(openxlsx)

# I/O
dat_file <- "output_data/out_val/table_val_guideScores_whole.txt"
deg_file <- "~/projects/GIN_2020/final/output/global_network/figures/GIN_global_network_GI_density_standard.xlsx"
qgi_file <- "~/data/qGI/20210107/qGI_20210107.txt"
ess_file <- "~/data/essentials/HAP1_core_essentials_predictive.txt"
dess_file <- "~/data/essentials/DepMap_essential_19Q2_60_percent.txt"

output_folder <- file.path(dirname(dat_file), "guidePosition")
if (!dir.exists(output_folder)) dir.create(output_folder)

# Read in
dat <- fread(dat_file, h = TRUE)
deg <- import_list(deg_file)
deg_all <- deg[["GI_library_density_all"]]
deg_all <- deg_all[,-5]
colnames(deg_all)[1] <- "gene_name"
qgi <- read.delim(qgi_file, h = TRUE)
queries <- sub("_[0-9][0-9][0-9]", "", colnames(qgi))
ess <- readLines(ess_file)
dess <- readLines(dess_file)

# Init list
hypo_list <- list()

# Loop through VBC score (guide efficiency prediction) thresholds
# and summarise tiling array LFC
for (i in c(0.90, 0.95)) {
    # Define vbc score threshold and summarise gene info
    quant <- quantile(dat$vbc_score, i, na.rm = TRUE)
    dat_vbc <- filter(dat, vbc_score >= quant)
    hypo_genes <- dat_vbc %>%
        mutate(position = ifelse(distance_tss_0_stop_1 <= 0.5, "tss", "end")) %>%
        mutate(position = factor(position, levels = c("tss", "end"))) %>%
        group_by(gene_name, position, .drop = FALSE) %>%
        summarise(n = n(), mean_LFC = mean(logFC, na.rm = TRUE), .groups = "keep") %>%
        pivot_wider(names_from = position, values_from = c(n, mean_LFC)) %>%
        drop_na() %>%
        mutate(deltaLFC = mean_LFC_tss - mean_LFC_end) %>%
        filter(deltaLFC < -2) %>%
        filter(n_tss + n_end >= 4) %>%
        left_join(deg_all, by = "gene_name") %>%
        mutate(depmap60_ess = ifelse(gene_name %in% dess, TRUE, FALSE)) %>%
        mutate(HAP1_ess = ifelse(gene_name %in% ess, TRUE, FALSE)) %>%
        mutate(screened = ifelse(gene_name %in% queries, TRUE, FALSE)) %>%
        as.data.frame()

    hypo_guides <- dat_vbc %>%
        filter(gene_name %in% hypo_genes$gene_name) %>%
        select(gene_name, chrom, start, stop, strand, guide, distance_tss_0_stop_1,
            guide_library, logFC, guide_strand, gene_strand, vbc_score, escapeNMD,
            provean_score, guide_target, Interpro_Description) %>%
        group_by(gene_name) %>%
        arrange(distance_tss_0_stop_1, .by_group = TRUE) %>%
        as.data.frame()

    # Populate list
    list_name1 <- paste0("vbc", i*100, "_", quant*100, "_genes")
    list_name2 <- paste0("vbc", i*100, "_", quant*100, "_guides")
    hypo_list[[list_name1]] <- hypo_genes
    hypo_list[[list_name2]] <- hypo_guides
}

# Write out
output_file <- file.path(output_folder, "tiling_vbc_position_deltaLFC2_predHypomorph.xlsx")
write.xlsx(hypo_list, output_file)
