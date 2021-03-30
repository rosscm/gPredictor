# Script to summarise guide feature scores across gene targeting positions
library(data.table)
library(dplyr)
library(ggplot2)
library(openxlsx)

# Read in guide score data
dat_file <- "output_data/out_val/table_val_guideScores_whole.txt"
dat <- fread(dat_file)
dat <- na.omit(dat)

output_folder <- file.path(dirname(dat_file), "guidePosition")
if (!dir.exists(output_folder)) dir.create(output_folder)

# Set distance threshold to define N/C terminus guides
quant1 <- quantile(dat$distance_tss_0_stop_1, 0.25)
quant3 <- quantile(dat$distance_tss_0_stop_1, 0.75)
dat$gposition <- "middle"
dat[which(dat$distance_tss_0_stop_1 <= quant1), "gposition"] <- "tss"
dat[which(dat$distance_tss_0_stop_1 >= quant3), "gposition"] <- "end"

# Prepare data for random forest classification
pred <- dat[,c("logFC", "Exon_Length", "vbc_score",
  "doench_score", "provean_score", "disorder_score",
  "guideEfficiency_pred", "predicted_in_frame", "predicted_oof",
  "sequence_score", "escapeNMD", "distance_tss_0_stop_1")]
pred <- na.omit(pred)

# Summarise scores for N/C terminus guide groups
dat_sum <- dat %>%
  group_by(gposition) %>%
  summarise(n_guides = n(),
            logFC = mean(logFC, na.rm = TRUE),
            Exon_Length = mean(Exon_Length, na.rm = TRUE),
            vbc_score = mean(vbc_score, na.rm = TRUE),
            #offtarget_score = mean(offtarget_score, na.rm = TRUE),
            doench_score = mean(doench_score, na.rm = TRUE),
            #CRISPRO_score = mean(CRISPRO_score, na.rm = TRUE),
            #oof_score = mean(oof_score, na.rm = TRUE),
            #num_transcripts = mean(num_transcripts, na.rm = TRUE),
            #targeted_transcripts = mean(targeted_transcripts_frac, na.rm = TRUE),
            provean_score = mean(provean_score, na.rm = TRUE),
            disorder_score = mean(disorder_score, na.rm = TRUE),
            #Frameshift_frequency = mean(Frameshift.frequency, na.rm = TRUE),
            #Expected_indel_length = mean(Expected.indel.length, na.rm = TRUE),
            #Precision = mean(Precision, na.rm = TRUE),
            guideEfficiency_pred = mean(guideEfficiency_pred, na.rm = TRUE),
            predicted_in_frame = mean(predicted_in_frame, na.rm = TRUE),
            predicted_oof = mean(predicted_oof, na.rm = TRUE),
            sequence_score = mean(sequence_score, na.rm = TRUE),
            #domain_target_perc = length(which(guide_target == "protein_domain")) / n() * 100,
            #protein_structure_perc = length(which(guide_structure == "protein_structure")) / n() * 100,
            escapeNMD_perc =  length(which(escapeNMD == TRUE)) / n() * 100
            #Exon_Multiple_of_3 = length(which(Exon_Multiple_of_3 == TRUE)) / n() * 100
          ) %>% as.data.frame

# Write out
table_out1 <- file.path(output_folder, "guide_feature_mean_NCterminus.xlsx")
write.xlsx(dat_sum, table_out1, row.names = FALSE)

# Plot
dat_melt <- reshape2::melt(dat_sum)
dat_melt$gposition <- factor(dat_melt$gposition, levels = c("tss", "middle", "end"))
p1 <- ggplot(dat_melt, aes(x = gposition, y = value)) +
      facet_wrap(. ~ variable, scales = "free") +
      geom_bar(stat = "identity") +
      labs(x = NULL, y = NULL, title = "Guide feature distribution across guide target positions",
          subtitle = sprintf("Tiling library (%s guides, %s genes)", nrow(dat), length(unique(dat$gene_name)))) +
      theme_linedraw() +
      theme(panel.grid = element_blank())

# Draw out
plot_out1 <- file.path(output_folder, "guide_feature_mean_NCterminus.pdf")
ggsave(plot_out1, p1, width = 7, height = 6)

# Positive vs. escapeNMD logFC
dat_sum2 <- dat %>%
  group_by(gposition, escapeNMD) %>%
  summarise(n_guides = n(),
            logFC_mean = mean(logFC, na.rm = TRUE),
            logFC_sd = sd(logFC, na.rm = TRUE)
          ) %>% as.data.frame

dat_sum2$gposition <- factor(dat_sum2$gposition, levels = c("tss", "middle", "end"))

# Plot
p2 <- ggplot(dat_sum2, aes(x = gposition, y = logFC_mean, fill = escapeNMD)) +
        geom_bar(stat = "identity", position = position_dodge()) +
        geom_text(aes(label = n_guides), position = position_dodge(width = 0.75), vjust = 1,size = 3) +
        labs(x = NULL, y = "Average dropout") +
        scale_fill_brewer(palette = "Dark2") +
        theme_linedraw() +
        theme(panel.grid = element_blank())

# Draw out
plot_out2 <- file.path(output_folder, "guide_feature_mean_NCterminus_logFC_NMD.pdf")
ggsave(plot_out2, p2, width = 4, height = 3)
