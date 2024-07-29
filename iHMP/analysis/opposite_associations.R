library(ggplot2)
library(gridExtra)
library(lmerTest)

taxa_table <- read.csv('iHMP/analysis/data/metaphlan4_taxonomic_profiles.tsv', skip = 1, check.names = F, sep = '\t')

# Reorganize taxa table
colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                           (grepl('\\|s__', taxa_table$clade_name) &
                              !grepl('\\|t__', taxa_table$clade_name)),]
rownames(taxa_table) <- taxa_table$clade_name
taxa_table$clade_name <- NULL

prepare_metadata <- function(dataset_type) {
  metadata <- read.csv('iHMP/analysis/data/hmp2_metadata_2018-08-20.csv', check.names = F)
  
  # Reorganize metadata table
  for (participant_id in unique(metadata$`Participant ID`)) {
    bmis <- metadata$BMI[metadata$`Participant ID` == participant_id]
    bmis <- bmis[!is.na(bmis)]
    metadata$BMI[metadata$`Participant ID` == participant_id] <- ifelse(length(bmis) > 0, mean(bmis), NA)
    
    smoke_status <- metadata$`smoking status`[metadata$`Participant ID` == participant_id]
    smoke_status <- smoke_status[!is.na(smoke_status)]
    metadata$`smoking status`[metadata$`Participant ID` == participant_id] <- ifelse(length(smoke_status) > 0, smoke_status[1], NA)
  }
  
  if (dataset_type == 'taxa') {
    metadata <- metadata[metadata$data_type == 'metagenomics',]
  } else {
    metadata <- metadata[metadata$data_type == 'metabolomics',]
  }
  rownames(metadata) <- metadata$`External ID`
  metadata <- metadata[,colSums(metadata == '', na.rm = T) != nrow(metadata)]
  keep_cols <- c('External ID', 'Participant ID', 'week_num', 'site_name', 'Age at diagnosis',
                 'Education Level', 'Occupation', 'consent_age', 'diagnosis',
                 colnames(metadata)[c(52:83, 85:111)], 'race', 'sex', 'BMI', 'reads_filtered')
  metadata <- metadata[,keep_cols]
  metadata <- metadata[,colSums(!is.na(metadata)) != 0]
  return(metadata)
}
metadata <- prepare_metadata('taxa')
metadata$age <- metadata$consent_age + metadata$week_num / 52

abun_df <- data.frame(sample = colnames(taxa_table),
                      abun = unlist(taxa_table['k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Lachnospiraceae_unclassified|s__Eubacterium_rectale',]))

# Plot 1: pseudo-counts
joined_df <- full_join(abun_df, metadata, by = c('sample' = 'External ID'))
joined_df$abun <- ifelse(joined_df$abun == 0, 
                         min(joined_df$abun[joined_df$abun > 0], na.rm=T)/2, 
                         joined_df$abun)

plot1 <- ggplot(joined_df, aes(x = age, y = abun)) + 
  geom_point(size = 0.5) + 
  scale_y_continuous(transform = 'log', breaks = 10^seq(-3, 1)) + 
  geom_smooth(method = 'lm') + 
  theme_bw() + 
  xlab("Age") + 
  ylab("Eubacterium Rectale\nRelative Abundance (%)")

summary(lm(log(abun) ~ age, joined_df))

# Plot 2: non-zero abundance
joined_df <- full_join(abun_df, metadata, by = c('sample' = 'External ID'))
joined_df$abun <- ifelse(joined_df$abun == 0, 
                          NA, 
                         joined_df$abun)

plot2 <- ggplot(joined_df, aes(x = age, y = abun)) + 
  geom_point(size = 0.5) + 
  scale_y_continuous(transform = 'log', breaks = 10^seq(-3, 1)) + 
  geom_smooth(method = 'lm') + 
  theme_bw() + 
  xlab("Age") + 
  ylab("Eubacterium Rectale\nRelative Abundance (%)")

summary(lm(log(abun) ~ age, joined_df))

# Plot 3: binary
joined_df <- full_join(abun_df, metadata, by = c('sample' = 'External ID'))
joined_df$abun <- ifelse(joined_df$abun == 0, 
                         0, 
                         1)

model <- glm(abun ~ age, data = joined_df, family = binomial)

x_seq <- seq(min(joined_df$age, na.rm=T), max(joined_df$age, na.rm=T), length.out = 1000)
pred <- predict(model, newdata = data.frame(age = x_seq), type = "response", se.fit = TRUE)
pred_df <- data.frame(
  x = x_seq,
  y = pred$fit,
  ymin = pred$fit - 1.96 * pred$se.fit,
  ymax = pred$fit + 1.96 * pred$se.fit
)

plot3 <- ggplot(joined_df, aes(x = age, y = abun)) + 
  geom_line(data = pred_df, aes(x = x, y = y), color = "blue", linewidth = 1) +
  geom_ribbon(data = pred_df, aes(x = x, y = y, ymin = ymin, ymax = ymax), alpha = 0.2) +
  geom_jitter(height = 0.05, size = 0.5) + 
  theme_bw() + 
  xlab("Age") + 
  ylab("Eubacterium Rectale\nPrevalence")

summary(glm(abun ~ age, joined_df, family = binomial))

combined_plot <- grid.arrange(plot1, plot2, plot3, ncol=3)
ggsave(combined_plot, filename = 'Figures/paper_figures/opposite_associations.png', width = 8, height = 2.5)

