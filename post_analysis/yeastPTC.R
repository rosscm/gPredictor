# Script to aggregate yeast PTC tolerance data and map to human orthologs
# https://www.nature.com/articles/s41588-018-0087-y#Sec21
library(openxlsx)
library(data.table)
library(dplyr)

# I/O
ortho_file <- "~/data/yeast/ORTHOLOGY-ALLIANCE_yeastHuman090320.xlsx"
ptc_file <- "~/data/yeast/yeast_hypomorph_PTC_paper.xlsx"
output_folder <- "output_data/out_val/guidePosition"
ortho <- read.xlsx(ortho_file)
ptc <- read.xlsx(ptc_file)
ptc <- data.table(ptc)


# Clean yeast-human ortholog data
# Pick best score for yeast genes with multiple human gene annotations
ortho2 <- ortho %>%
  select(Yeast.Gene, Human.Gene, IsBestScore) %>%
  arrange(Yeast.Gene) %>%
  group_by(Yeast.Gene) %>%
  subset(IsBestScore == "Yes") %>%
  aggregate(Human.Gene ~ Yeast.Gene, ., function(x) paste(x, collapse = ", ")) %>%
  rename(gene_name = Yeast.Gene) %>%
  rename(human_gene_name = Human.Gene)

# Select informative columns from PTC data and join with ortho2
ptc_all <- ptc %>%
  select(gene_name, guide, CDS_length, dist_from_CDS_end,
      "PTC.tolerance.score.(combined)", "gene.PTC.tolerance.(combined)",
      dubious, viable_annotation, "evolvability.(PMID:.26627736)", "human.complements.null.(PMID:.25999509)") %>%
  mutate(gene_name = gsub(" ", "", gene_name)) %>%
  left_join(ortho2, by = "gene_name")
colnames(ptc_all) <- gsub("\\.|\\(|\\)|\\:", "_", colnames(ptc_all))

# Filter for yeast genes with human orthologs
ptc_sub1 <- ptc_all %>%
  filter(!is.na(human_gene_name)) %>%
  arrange(gene_name, dist_from_CDS_end)

# Filter for genes with positive tolerance score (can tolerate pre-mature stop codon),
ptc_sub2 <- ptc_sub1 %>%
  filter(gene_PTC_tolerance__combined_ > 1)

# Prepare output list
ptc_list <- list()
ptc_list[["all_data"]] <- ptc
ptc_list[["human_ortho"]] <- ptc_sub1
ptc_list[["human_ortho_genePTC1"]] <- ptc_sub2

# Write out
fname <- file.path(output_folder, "yeast_human_ptc_tolerance_guides.xlsx")
write.xlsx(ptc_list, file = fname, row.names = FALSE)
