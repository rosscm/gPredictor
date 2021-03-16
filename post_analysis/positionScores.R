# Script to summarise guide feature scores across gene targeting positions
library(dplyr)
library(ggplot2)

dat <- read.delim("output_data/out_vbc/table_val_guideScores_whole.txt", h = TRUE, as.is = TRUE)
distance90 <- quantile(dat$distance_tss_0_stop_1, 0.90, na.rm=TRUE)
dat$tss_end <- "middle"
dat[which(dat$distance_tss_0_stop_1 < 1-distance90), "tss_end"] <- "tss"
dat[which(dat$distance_tss_0_stop_1 > distance90), "tss_end"] <- "end"

dat_sum <- dat %>%
  group_by(tss_end, .drop=T) %>%
  summarise(n_guides = n(),
            logFC = mean(logFC, na.rm = TRUE),
            Exon_Length = mean(Exon_Length, na.rm = TRUE),
            vbc_score = mean(vbc_score, na.rm = TRUE),
            offtarget_score = mean(offtarget_score, na.rm = TRUE),
            doench_score = mean(doench_score, na.rm = TRUE),
            CRISPRO_score = mean(CRISPRO_score, na.rm = TRUE),
            oof_score = mean(oof_score, na.rm = TRUE),
            num_transcripts = mean(num_transcripts, na.rm = TRUE),
            targeted_transcripts = mean(targeted_transcripts_frac, na.rm = TRUE),
            provean_score = mean(provean_score, na.rm = TRUE),
            disorder_score = mean(disorder_score, na.rm = TRUE),
            Frameshift_frequency = mean(Frameshift.frequency, na.rm = TRUE),
            Expected_indel_length = mean(Expected.indel.length, na.rm = TRUE),
            Precision = mean(Precision, na.rm = TRUE),
            guideEfficiency_pred = mean(guideEfficiency_pred, na.rm = TRUE),
            predicted_in_frame = mean(predicted_in_frame, na.rm = TRUE),
            predicted_oof = mean(predicted_oof, na.rm = TRUE),
            sequence_score = mean(sequence_score, na.rm = TRUE),
            domain_target_perc = length(which(guide_target == "protein_domain")) / n() * 100,
            protein_structure_perc = length(which(guide_structure == "protein_structure")) / n() * 100,
            SecStruct_perc = length(which(SecStruct != "None")) / n() * 100,
            escapeNMD_perc =  length(which(escapeNMD == TRUE)) / n() * 100,
            Exon_Multiple_of_3 = length(which(Exon_Multiple_of_3 == TRUE)) / n() * 100
          ) %>% as.data.frame

dat2 <- reshape2::melt(dat_sum)
dat2$tss_end <- factor(dat2$tss_end, levels = c("tss", "middle", "end"))

ggplot(dat2, aes(x = tss_end, y = value)) +
  facet_wrap(. ~ variable, scales = "free") +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = NULL, title = "Guide feature distribution across guide target positions",
      subtitle = sprintf("Tiling library (%s guides, %s genes)", nrow(dat), length(unique(dat$gene_name)))) +
  theme_linedraw() +
  theme(panel.grid = element_blank())
