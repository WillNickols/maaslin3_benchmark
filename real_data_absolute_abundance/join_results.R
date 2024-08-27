library(UpSetR)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggExtra)
library(gridExtra)
library(ggbeeswarm)

figures_folder <- paste0('Figures/paper_figures/')

all_results <- c(list.files('real_data_absolute_abundance/Barlow/analysis/results/', full.names = T),
                 list.files('real_data_absolute_abundance/infants/analysis/results/', full.names = T),
                 list.files('real_data_absolute_abundance/VieiraSilva/analysis/results/', full.names = T))
                 # list.files('real_data_absolute_abundance/nguyen_covid/analysis/results/', full.names = T))

problematic_categoricals <- c(DietKeto_DietKeto = 'Diet_Keto',
                              diagnosisCD_diagnosisCD = 'diagnosis_CD',
                              diagnosisPSC_diagnosisPSC = 'diagnosis_PSC',
                              `diagnosisPSC-CD_diagnosisPSC-CD` = 'diagnosis_PSC-CD',
                              `diagnosisPSC-UC_diagnosisPSC-UC` = 'diagnosis_PSC-UC',
                              genderM_genderM = 'gender_M',
                              covid_19_severitySevere_covid_19_severitySevere = 'covid_19_severity_Severe',
                              ethnicityhispanic_ethnicityhispanic = 'ethnicity_hispanic',
                              raceAsian_raceAsian = 'race_Asian',
                              raceBlack_raceBlack = 'race_Black',
                              raceOther_raceOther = 'race_Other',
                              remdesiviryes_remdesiviryes = 'remdesivir_yes',
                              corticosteroidsyes_corticosteroidsyes = 'corticosteroids_yes',
                              sexfemale_sexfemale = 'sex_female',
                              antibioticsyes_antibioticsyes = 'antibiotics_yes')

growing_df <- data.frame()
for (result in all_results) {
  fit_out_joint <- read.csv(result, sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
  
  if (grepl('Maaslin2', result)) {
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$tool <- 'Maaslin2'
  }
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint <- fit_out_joint[fit_out_joint$association == 'abundance',]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_joint
    fit_out_joint$qval <- fit_out_joint$qval_joint
    fit_out_joint$tool <- gsub('\\.tsv', '', gsub('.*_', '', result))
  }
  
  if (grepl('ANCOMBC', result)) {
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error) | fit_out_joint$error == 'sensitivity failed',] # Allow sensitivity failed because we just care about the coefficients
    fit_out_joint <- fit_out_joint[!is.infinite(fit_out_joint$effect_size),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, "_", fit_out_joint$metadata)
    fit_out_joint$metadata_value[fit_out_joint$metadata_value %in% names(problematic_categoricals)] <- 
      problematic_categoricals[fit_out_joint$metadata_value[fit_out_joint$metadata_value %in% names(problematic_categoricals)]]
    fit_out_joint$tool <- 'ANCOMBC'
    fit_out_joint$coef <- fit_out_joint$effect_size / log(2)
  }
  
  if (grepl('ALDEx2', result)) {
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, "_", fit_out_joint$metadata)
    fit_out_joint$metadata_value[fit_out_joint$metadata_value %in% names(problematic_categoricals)] <- 
      problematic_categoricals[fit_out_joint$metadata_value[fit_out_joint$metadata_value %in% names(problematic_categoricals)]]
    fit_out_joint$tool <- 'ALDEx2'
    fit_out_joint$coef <- fit_out_joint$effect_size
    fit_out_joint <- fit_out_joint[!grepl('mouse', fit_out_joint$metadata_value),]
  }
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool')]
  fit_out_joint$study <- gsub("real_data_absolute_abundance\\/|\\/analysis.*", "", result)
  
  growing_df <- rbind(growing_df, fit_out_joint)
}

names_excluded <- c('read_depth_read_depth', 
                    'read_count_read_count', 
                    'remdesivir_yes', 
                    'corticosteroids_yes',
                    'race_Asian',
                    'race_Black',
                    'race_Other',
                    'ethnicity_hispanic',
                    'sex_female', 
                    'antibiotics_yes',
                    'gender_M')

growing_df <- growing_df[!growing_df$metadata_value %in% names_excluded,]
growing_df <- growing_df %>%
  mutate(metadata_value = case_when(metadata_value == 'age_age' ~ "Age",
                                    metadata_value == 'bmi_bmi' ~ "BMI",
                                    metadata_value == 'day_day' ~ "Day",
                                    metadata_value == 'Day_Day' ~ "Day",
                                    metadata_value == 'diagnosis_CD' ~ "CD",
                                    metadata_value == 'diagnosis_PSC' ~ "PSC",
                                    metadata_value == 'diagnosis_PSC-CD' ~ "PSC-CD",
                                    metadata_value == 'diagnosis_PSC-UC' ~ "PSC-UC",
                                    metadata_value == 'Diet_Keto' ~ "Keto Diet",
                                    metadata_value == 'gender_M' ~ "Male",
                                    metadata_value == 'BMI_BMI' ~ "BMI",
                                    metadata_value == 'cci_cci' ~ "CCI",
                                    metadata_value == 'covid_19_severity_Severe' ~ "COVID-19 Severity",))

custom_colors <- c("Infant gut: Day" = "deeppink1",
                   "Mouse gut: Day" = "#FF7F00",
                   "Mouse gut: Keto Diet" = "green4",
                   "IBD/PSC gut: Age" = "dodgerblue2", 
                   "IBD/PSC gut: BMI" = "#6A3D9A", 
                   # "IBD/PSC gut: Male" = "#E31A1C", 
                   "IBD/PSC gut: CD" = "palegreen2", 
                   "IBD/PSC gut: PSC" = "gray70", 
                   "IBD/PSC gut: PSC-CD" = "gold1", 
                   "IBD/PSC gut: PSC-UC" = "skyblue2",
                   "COVID-19: Age" = "blue1",
                   "COVID-19: BMI" = "#FB9A99",
                   "COVID-19: CCI" = "#CAB2D6",
                   "COVID-19: COVID-19 Severity" = "khaki2")

plot_list <- list()
for (study in unique(growing_df$study)) {
  df_tmp_1 <- growing_df[growing_df$tool == 'Maaslin3Inferred' & growing_df$study == study,]
  df_tmp_2 <- growing_df[growing_df$tool == 'Maaslin3CompAdjust' & growing_df$study == study,]
  
  tmp_join <- full_join(df_tmp_1, df_tmp_2, by = c('feature', 'metadata_value'))
  colnames(tmp_join) <- gsub('\\.x', '_Absolute', colnames(tmp_join))
  colnames(tmp_join) <- gsub('\\.y', '_Other', colnames(tmp_join))
  tmp_join <- tmp_join[!is.na(tmp_join$coef_Absolute) & !is.na(tmp_join$coef_Other),]
  tmp_join$coef_Absolute[tmp_join$coef_Absolute > median(tmp_join$coef_Absolute) + 
                           (quantile(tmp_join$coef_Absolute, 0.98) - median(tmp_join$coef_Absolute)) * 2] <- Inf
  tmp_join$coef_Absolute[tmp_join$coef_Absolute < median(tmp_join$coef_Absolute) + 
                           (quantile(tmp_join$coef_Absolute, 0.02) - median(tmp_join$coef_Absolute)) * 2] <- -Inf
  tmp_join$coef_Other[tmp_join$coef_Other > median(tmp_join$coef_Other) + 
                           (quantile(tmp_join$coef_Other, 0.98) - median(tmp_join$coef_Other)) * 2] <- Inf
  tmp_join$coef_Other[tmp_join$coef_Other < median(tmp_join$coef_Other) + 
                           (quantile(tmp_join$coef_Other, 0.02) - median(tmp_join$coef_Other)) * 2] <- -Inf
  
  study_name <- case_when(study == 'Barlow' ~ "\nMouse gut",
                          study == 'infants' ~ "\nInfant gut",
                          study == 'VieiraSilva' ~ "\nIBD/PSC gut",
                          study == 'nguyen_covid' ~ "\nCOVID-19")
  
  tmp_join$meta_in_study <- paste0(gsub('\\\n', '', study_name), ': ', tmp_join$metadata_value)
  
  p <- ggplot(data = tmp_join, aes(x = coef_Absolute, y = coef_Other, color = meta_in_study)) + 
    geom_point(aes(colour=meta_in_study)) + 
    geom_point(shape = 21, colour = "black", aes(fill = meta_in_study)) +
    geom_abline(slope = 1, intercept = 0) + 
    geom_vline(xintercept = 0, linewidth = 0.2) + 
    geom_hline(yintercept = 0, linewidth = 0.2) + 
    theme_bw()+ 
    xlab('Absolute coefficient') + 
    ylab('Relative coefficient') + 
    theme(legend.position = 'none') + 
    labs(color = '') +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    theme(text = element_text(size = 16)) + 
    ggtitle(study_name)
  
  plot_list[[study]] <- ggMarginal(p, type = "histogram", margins = "both", groupFill = T)
}

plot_out <- grid.arrange(plot_list[[2]], plot_list[[1]], plot_list[[3]], ncol = 3)
ggsave(paste0(figures_folder, 'abs_vs_rel_real_data.png'),
       plot = plot_out, width = 12, height = 5)

growing_df_2 <- data.frame()
for (tool in unique(growing_df$tool)) {
  if (tool == 'Maaslin3Inferred') {
    next
  }
  df_tmp_1 <- growing_df[growing_df$tool == 'Maaslin3Inferred',]
  df_tmp_2 <- growing_df[growing_df$tool == tool,]
  
  tmp_join <- full_join(df_tmp_1, df_tmp_2, by = c('feature', 'metadata_value', 'study'))
  colnames(tmp_join) <- gsub('\\.x', '_Absolute', colnames(tmp_join))
  colnames(tmp_join) <- gsub('\\.y', '_Other', colnames(tmp_join))
  tmp_join <- tmp_join[!is.na(tmp_join$coef_Absolute) & !is.na(tmp_join$coef_Other),]
  tmp_join <- tmp_join %>% 
    group_by(study) %>%
    mutate(coef_Absolute = ifelse(coef_Absolute > median(coef_Absolute) + (quantile(coef_Absolute, 0.98) - median(coef_Absolute)) * 2 | 
                                    coef_Absolute < median(coef_Absolute) + (quantile(coef_Absolute, 0.02) - median(coef_Absolute)) * 2, NA, coef_Absolute),
           coef_Other = ifelse(coef_Other > median(coef_Other) + (quantile(coef_Other, 0.98) - median(coef_Other)) * 2 | 
                                 coef_Other < median(coef_Other) + (quantile(coef_Other, 0.02) - median(coef_Other)) * 2, NA, coef_Other))

  lm_df <- tmp_join %>%
    dplyr::group_by(metadata_value, study) %>%
    dplyr::summarize(coef(lm(coef_Other ~ coef_Absolute))[-1]) # Get only the slope
  colnames(lm_df) <- c('metadata_value', 'study', 'lm')
  
  cor_df <- tmp_join %>%
    dplyr::group_by(metadata_value, study) %>%
    dplyr::summarize(cor(coef_Absolute,
                  coef_Other,
                  use = 'pairwise.complete.obs',
                  method = 'spearman'))
  colnames(cor_df) <- c('metadata_value', 'study', 'cor')
  
  joined_df <- full_join(lm_df, cor_df, by = c('metadata_value', 'study'))
  joined_df$tool <- tool
  growing_df_2 <- rbind(growing_df_2, joined_df)
}

growing_df_2 <- melt(growing_df_2, id.vars = c('metadata_value', 'study', 'tool'))

growing_df_2 <- growing_df_2 %>%
  mutate(tool = case_when(tool == 'ALDEx2' ~ 'ALDEx2',
                         tool == 'ANCOMBC' ~ 'ANCOM-BC2',
                         tool == 'Maaslin2' ~ 'MaAsLin 2',
                         tool == 'Maaslin3CompAdjust' ~ 'MaAsLin 3'),
         variable = case_when(variable == 'lm' ~ 'Slope (Fit vs. Absolute)',
                              variable == 'cor' ~ 'Abundance coefficient correlation'),
         study = case_when(study == 'Barlow' ~ "Mouse gut",
                           study == 'infants' ~ "Infant gut",
                           study == 'VieiraSilva' ~ "IBD/PSC gut",
                           study == 'nguyen_covid' ~ "COVID-19",),
         meta_in_study = factor(paste0(study, ': ', metadata_value), 
                                levels = names(custom_colors))
  )

plot_out_2 <- ggplot(growing_df_2, aes(x = tool, y = value, fill = meta_in_study)) + 
  geom_beeswarm(aes(shape = variable), size = 3, cex = 3.3, shape = 21) + 
  facet_wrap(. ~ variable) + 
  scale_y_continuous(breaks = seq(0, 1, 0.2)) + 
  theme_bw() + 
  xlab(NULL) + 
  ylab(NULL) + 
  scale_fill_manual(values = custom_colors) +
  guides(shape = "none") + 
  labs(fill = 'Metadatum') + 
  theme(text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        strip.background = element_rect(fill = "gray95"),
        panel.grid.minor = element_blank())
ggsave(paste0(figures_folder, 'abs_vs_real_all_tools.png'),
       plot = plot_out_2, width = 10, height = 3.5)

# In-text numbers
growing_df_2 %>%
    filter(tool == "MaAsLin 3")

growing_df_2 %>%
    dplyr::filter(study != "Infant gut") %>%
    dplyr::group_by(tool, variable) %>%
    dplyr::summarize(mean(value))

taxa_table <- read.csv('real_data_absolute_abundance/infants/analysis/data/41586_2021_3241_MOESM4_ESM.csv', skip = 1, check.names = F)
taxa_table <- taxa_table[,grepl("NICU|OTU_ID", colnames(taxa_table))]
taxa_table <- taxa_table[rowSums(is.na(taxa_table)) == 0,]

rownames(taxa_table) <- taxa_table$OTU_ID
taxa_table$OTU_ID <- NULL

mean(taxa_table == 0)

# Read in data
taxa_table <- read.csv('real_data_absolute_abundance/Barlow/analysis/data/Absolute_Abundance_Table.csv', check.names = F, sep = ',')
rownames(taxa_table) <- taxa_table[,1]
taxa_table[,1] <- NULL

taxa_table <- taxa_table[,!colnames(taxa_table) %in% c('Diet', 'Site', 'Day', 'mouse', 'Cage')]
colnames(taxa_table) <- make.names(colnames(taxa_table))
mean(taxa_table == 0)

taxa_table <- read.csv('real_data_absolute_abundance/VieiraSilva/analysis/data/QMP.matrix.tsv', check.names = F, sep = '\t')
rownames(taxa_table) <- taxa_table[,1]
taxa_table[,1] <- NULL
taxa_table <- t(apply(taxa_table, 1, function(x) {x / sum(x)})) # Convert to relative abundance
mean(taxa_table == 0)



