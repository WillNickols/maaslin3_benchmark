library(UpSetR)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)

workingDirectory <- "~/Documents/GitHub/maaslin3_benchmark"
figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')

taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', sep = '\t', skip = 1)

colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                           (grepl('\\|s__', taxa_table$clade_name) &
                              !grepl('\\|t__', taxa_table$clade_name)),]
rownames(taxa_table) <- taxa_table$clade_name
taxa_table$clade_name <- NULL
keep_taxa <- gsub(".*s__", "", rownames(taxa_table)[rowMeans(taxa_table > 0.1) > 0.01])

all_results <- list.files('results/')

growing_df <- data.frame()
for (result in all_results[grepl('^ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
  
  if (grepl('Maaslin2', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$tool <- 'Maaslin2'
  }
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_joint
    fit_out_joint$qval <- fit_out_joint$qval_joint
    if (grepl('Maaslin3ItAug', result)) {
      fit_out_joint$tool <- 'Maaslin3ItAug'
    } else {
      fit_out_joint$tool <- 'Maaslin3'
    }
  }
  
  if (grepl('ANCOMBC', result)) {
    fit_out_joint <- fit_out_joint[grepl('diagnosis|dysbiosis_state', fit_out_joint$metadata),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint <- fit_out_joint[!is.infinite(fit_out_joint$effect_size),]
    fit_out_joint$metadata_value <- gsub('dysbiosis_state|diagnosis', '', fit_out_joint$metadata)
    fit_out_joint$metadata_value <- ifelse(fit_out_joint$metadata_value %in% c('UC', 'CD'), 
                                           paste0('diagnosis_', fit_out_joint$metadata_value),
                                           paste0('dysbiosis_state_', fit_out_joint$metadata_value))
    fit_out_joint$tool <- 'ANCOMBC'
    fit_out_joint$coef <- fit_out_joint$effect_size
  }
  
  if (grepl('ALDEx2', result)) {
    fit_out_joint <- fit_out_joint[grepl('diagnosis|dysbiosis_state', fit_out_joint$metadata),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- gsub('dysbiosis_state|diagnosis', '', fit_out_joint$metadata)
    fit_out_joint$metadata_value <- ifelse(fit_out_joint$metadata_value %in% c('UC', 'CD'), 
                                           paste0('diagnosis_', fit_out_joint$metadata_value),
                                           paste0('dysbiosis_state_', fit_out_joint$metadata_value))
    fit_out_joint$tool <- 'ALDEx2'
    fit_out_joint$coef <- fit_out_joint$effect_size
  }
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool')]
  
  growing_df <- rbind(growing_df, fit_out_joint)
}

growing_df <- growing_df[growing_df$qval < 0.1 & growing_df$feature %in% keep_taxa,]
growing_df$association <- paste0(growing_df$feature, '-', growing_df$metadata_value)

listInput <- list(`Maaslin 2` = growing_df$association[growing_df$tool == 'Maaslin2'], 
     `Maaslin 3 Base` = growing_df$association[growing_df$tool == 'Maaslin3'], 
     `Maaslin 3` = growing_df$association[growing_df$tool == 'Maaslin3ItAug'], 
     `ANCOM-BC2` = growing_df$association[growing_df$tool == 'ANCOMBC'],
     ALDEx2 = growing_df$association[growing_df$tool == 'ALDEx2'])

plot_out <- upset(fromList(listInput), order.by = "freq", text.scale = 2.5)

png(file=paste0(figures_folder, 'fig_13a.png'), width = 14, height = 6, res = 1000, units = 'in')
plot_out
dev.off()

# Scatter plot
growing_df <- data.frame()
for (result in all_results[grepl('Maaslin3', all_results) & grepl('^ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_single
    fit_out_joint$qval <- fit_out_joint$qval_single
    if (grepl('Maaslin3ItAug', result)) {
      fit_out_joint$tool <- 'Maaslin3ItAug'
    } else {
      fit_out_joint$tool <- 'Maaslin3'
    }
  }
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool', 'association')]
  
  growing_df <- rbind(growing_df, fit_out_joint)
}

growing_df <- growing_df[growing_df$feature %in% keep_taxa,]
growing_df <- growing_df[growing_df$feature %in% growing_df$feature[growing_df$tool == 'Maaslin3ItAug'],]
growing_df$association_name <- paste0(growing_df$feature, '-', growing_df$metadata_value, '-', growing_df$association)

for_scatter <- growing_df[growing_df$association == 'abundance', c("association_name", "coef", "qval", "tool")]
split1 <- for_scatter[for_scatter$tool == 'Maaslin3',]
split2 <- for_scatter[for_scatter$tool == 'Maaslin3ItAug',]
colnames(split1) <- c('association_name', 'Maaslin3', 'qval_Maaslin3', 'tool')
colnames(split2) <- c('association_name', 'Maaslin3ItAug', 'qval_Maaslin3ItAug', 'tool')
split1$tool <- NULL
split2$tool <- NULL
joined_df <- left_join(split2, split1, by=c('association_name'))
joined_df$color <- ifelse(joined_df$qval_Maaslin3ItAug < 0.1,
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Both Significant', 
                                 'Maaslin 3 Only'),
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Maaslin 3 Base Only', 'FF'))
joined_df <- joined_df[!is.na(joined_df$color) & joined_df$color != 'FF',]

plot_out <- ggplot(joined_df, aes(x = Maaslin3ItAug, 
                      y = Maaslin3, color = color)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  geom_abline(slope=1, intercept=0, linetype = 'dashed') + 
  labs(color = '', x = 'Maaslin 3 Coefficient',
       y = 'Maaslin 3 Base Coefficient') + 
  scale_color_brewer(palette="Set1")
ggsave(paste0(figures_folder, 'fig_13c.png'),
       plot = plot_out, width = 9, height = 5)

for_scatter <- growing_df[growing_df$association == 'prevalence', c("association_name", "coef", "qval", "tool")]
split1 <- for_scatter[for_scatter$tool == 'Maaslin3',]
split2 <- for_scatter[for_scatter$tool == 'Maaslin3ItAug',]
colnames(split1) <- c('association_name', 'Maaslin3', 'qval_Maaslin3', 'tool')
colnames(split2) <- c('association_name', 'Maaslin3ItAug', 'qval_Maaslin3ItAug', 'tool')
split1$tool <- NULL
split2$tool <- NULL
joined_df <- left_join(split2, split1, by=c('association_name'))
joined_df$color <- ifelse(joined_df$qval_Maaslin3ItAug < 0.1,
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Both Significant', 
                                 'Maaslin 3 Only'),
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Maaslin 3 Base Only', 'FF'))
joined_df <- joined_df[!is.na(joined_df$color) & joined_df$color != 'FF',]

plot_out <- ggplot(joined_df, aes(x = Maaslin3ItAug, 
                                  y = Maaslin3, color = color)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  geom_abline(slope=1, intercept=0, linetype = 'dashed') + 
  labs(color = '', x = 'Maaslin 3 Coefficient',
       y = 'Maaslin 3 Base Coefficient') + 
  scale_color_brewer(palette="Set1")
ggsave(paste0(figures_folder, 'fig_13b.png'),
       plot = plot_out, width = 9, height = 5)

# Volcano plot
fit_out_joint <- read.csv(paste0('results/ibd_associations_Maaslin3ItAug.tsv'), sep = '\t')
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
fit_out_joint <- fit_out_joint[fit_out_joint$feature %in% keep_taxa,]

fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
fit_out_joint$pval <- fit_out_joint$pval_single
fit_out_joint$qval <- fit_out_joint$qval_single

# In-text number
mean(abs(fit_out_joint$coef[fit_out_joint$qval_single < 0.1]) < 10) - 
  mean(abs(fit_out_joint$coef[fit_out_joint$qval_single < 0.1]) < 1)

mean(abs(fit_out_joint$coef[fit_out_joint$qval_single < 0.1]) < 5) - 
  mean(abs(fit_out_joint$coef[fit_out_joint$qval_single < 0.1]) < 2)

significance_threshold <- 0.1
fit_out_joint$color <- ifelse(fit_out_joint$qval < significance_threshold, 
                              ifelse(fit_out_joint$association == 'prevalence', 'Significant prevalence', 'Significant abundance'),
                              ifelse(fit_out_joint$association == 'prevalence', 'Insignificant prevalence', 'Insignificant abundance'))
volcano_plot <- ggplot(fit_out_joint[sample(1:nrow(fit_out_joint)),], aes(x = coef, y = -log10(qval))) +
  geom_point(aes(color = as.factor(color)), size = 3, alpha=0.5) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
  labs(x = "Effect Size", y = "-log10(q-value)") +
  theme_bw() + 
  theme(text = element_text(size = 35),
        legend.position = 'right',
        legend.direction = 'vertical') +
  labs(color = '') + 
  scale_color_manual(values=c("#ff9496", "#9bd0e6", "#f71c1a", "#007ab9"))

ggsave(paste0(figures_folder, 'fig_13d.png'),
       plot = volcano_plot, width = 15, height = 6)

top_hits <- fit_out_joint[fit_out_joint$association == 'abundance' & fit_out_joint$feature != 'UNCLASSIFIED',]
top_hits <- rbind(top_hits[top_hits$coef > 0,][order(top_hits[top_hits$coef > 0,]$qval_single)[1:8],],
                  top_hits[top_hits$coef < 0,][order(top_hits[top_hits$coef < 0,]$qval_single)[1:8],])
top_hits <- rbind(top_hits[order(-top_hits$coef)[1:8],],
                  top_hits[order(top_hits$coef)[8:1],])
top_hits$feature <- gsub('_', ' ', top_hits$feature)
top_hits$feature <- paste0(top_hits$feature, ' - ', str_to_title(gsub("_", " ", top_hits$metadata)), ": ", gsub("dysbiosis ", "", gsub("_", " ", top_hits$value)))
top_hits$feature <- factor(top_hits$feature, levels = top_hits$feature)

labs1 <- levels(top_hits$feature)
labs2 <- case_when(top_hits$qval_single < 0.001 ~ "***",
                   top_hits$qval_single < 0.01 ~ "**",
                   top_hits$qval_single < 0.1 ~ "*")

plot_out <- ggplot(top_hits, aes(x = as.numeric(feature), y = coef)) + 
  geom_point(size = 3) + 
  coord_flip() + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  labs(x = '', y = 'Abundance coefficient') + 
  scale_x_continuous(breaks = 1:length(labs1),
                     labels = labs1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labs2),
                                         labels = labs2))

ggsave(paste0(figures_folder, 'fig_14a.png'),
       plot = plot_out, width = 14, height = 6)

top_hits <- fit_out_joint[fit_out_joint$association == 'prevalence' & fit_out_joint$feature != 'UNCLASSIFIED',]
top_hits <- rbind(top_hits[top_hits$coef > 0,][order(top_hits[top_hits$coef > 0,]$qval_single)[1:8],],
                  top_hits[top_hits$coef < 0,][order(top_hits[top_hits$coef < 0,]$qval_single)[1:8],])
top_hits <- rbind(top_hits[order(-top_hits$coef)[1:8],],
                  top_hits[order(top_hits$coef)[8:1],])
top_hits$feature <- gsub('_', ' ', top_hits$feature)
top_hits$feature <- paste0(top_hits$feature, ' - ', str_to_title(gsub("_", " ", top_hits$metadata)), ": ", gsub("dysbiosis ", "", gsub("_", " ", top_hits$value)))
top_hits$feature <- factor(top_hits$feature, levels = top_hits$feature)

labs1 <- levels(top_hits$feature)
labs2 <- case_when(top_hits$qval_single < 0.001 ~ "***",
                   top_hits$qval_single < 0.01 ~ "**",
                   top_hits$qval_single < 0.1 ~ "*")

plot_out <- ggplot(top_hits, aes(x = as.numeric(feature), y = coef)) + 
  geom_point(size = 3) + 
  coord_flip() + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  labs(x = '', y = 'Prevalence coefficient') + 
  scale_x_continuous(breaks = 1:length(labs1),
                     labels = labs1,
                     sec.axis = sec_axis(~.,
                                         breaks = 1:length(labs2),
                                         labels = labs2))

ggsave(paste0(figures_folder, 'fig_14b.png'),
       plot = plot_out, width = 14, height = 6)

############
# OMP/GOMP #
############

growing_df <- data.frame()
for (result in all_results[grepl('^omp|^gomp', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature) & fit_out_joint$metadata != "participant_id",]
  fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
  
  fit_out_joint <- fit_out_joint[grepl('food_group', fit_out_joint$metadata),]
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
  fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
  fit_out_joint$pval <- fit_out_joint$pval_single
  fit_out_joint$qval <- fit_out_joint$qval_single
  fit_out_joint$tool <- 'Maaslin3'
  fit_out_joint$method <- gsub("_.*", "", result)
  fit_out_joint$food <- gsub('\\.tsv$', '', gsub(".*associations_", "", result))
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'association', 'method', 'food')]
  
  growing_df <- rbind(growing_df, fit_out_joint[!is.na(fit_out_joint$qval),])
}

growing_df$qval <- p.adjust(growing_df$pval, method = 'BH')

signif_df <- growing_df[growing_df$qval < 0.1 & growing_df$feature %in% keep_taxa,]

results_table <- table(signif_df$feature, signif_df$method)
results_table <- results_table[rowSums(results_table) > 1,]
melted_df <- melt(results_table)
melted_df$Var1 <- gsub("_", " ", melted_df$Var1)
melted_df$Var1 <- factor(melted_df$Var1, levels = gsub("_", " ", rownames(results_table)[order(-rowSums(results_table))]))
melted_df$Var2 <- ifelse(melted_df$Var2 == 'gomp', "GOMP", "OMP")
plot_out <- ggplot(melted_df, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size = 15),
        legend.position = 'bottom') + 
  labs(fill = 'Association type') + 
  ylab('Significant calls across all diet groups') + 
  xlab('') + 
  scale_fill_brewer(palette="Accent")
ggsave(paste0(figures_folder, 'fig_15a.png'),
       plot = plot_out, width = 15, height = 8)

results_table <- table(signif_df$food, signif_df$method)
results_table <- results_table[rowSums(results_table) > 1,]
melted_df <- melt(results_table)
melted_df$Var1 <- factor(str_to_title(gsub('_', ' ', melted_df$Var1)), levels = str_to_title(gsub('_', ' ', rownames(results_table)[order(-rowSums(results_table))])))
melted_df$Var2 <- ifelse(melted_df$Var2 == 'gomp', "GOMP", "OMP")
plot_out <- ggplot(melted_df, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(position="stack", stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  labs(fill = 'Association type') + 
  ylab('Significant calls across\nall associations') + 
  xlab('') + 
  theme(text = element_text(size = 15)) + 
  scale_fill_brewer(palette="Accent")
ggsave(paste0(figures_folder, 'fig_15b.png'),
       plot = plot_out, width = 15, height = 4)

# results_table <- table(signif_df$metadata_value, signif_df$method)
# results_table <- results_table[rowSums(results_table) > 1,]
# melted_df <- melt(results_table)
# melted_df$Var1 <- gsub('food_group_', '', as.character(melted_df$Var1))
# melted_df$Var1 <- ifelse(melted_df$Var1 == 'food_group', 'Food group', melted_df$Var1)
# melted_df$Var1 <- factor(melted_df$Var1, levels = c('Food group', 
#                                                     'Within the past 4 to 7 days',
#                                                     'Within the past 2 to 3 days',
#                                                     'Yesterday, 1 to 2 times',
#                                                     'Yesterday, 3 or more times'))
# ggplot(melted_df, aes(x = Var1, y = value, fill = Var2)) + 
#   geom_bar(position="stack", stat="identity") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
#   labs(fill = 'Predictor type') + 
#   ylab('Significant associations across all food') + 
#   xlab('') + 
#   theme(text = element_text(size = 15))
#   
# signif_df[signif_df$feature == 'Lachnospiraceae_bacterium',]
# 
# signif_df$color <- ifelse(signif_df$qval < significance_threshold, 
#                               ifelse(signif_df$association == 'prevalence', 'Significant prevalence', 'Significant abundance'),
#                               ifelse(signif_df$association == 'prevalence', 'Insignificant prevalence', 'Insignificant abundance'))
# volcano_plot <- ggplot(signif_df, aes(x = coef, y = -log10(qval))) +
#   geom_point(aes(color = as.factor(metadata_value)), size = 3, alpha=0.5) +
#   geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
#   labs(x = "Effect Size", y = "-log10(q-value)", color = 'Metadata') +
#   theme_bw() + 
#   theme(text = element_text(size = 25))
# 
# volcano_plot

growing_df <- data.frame()
for (result in all_results[grepl('Maaslin3|Maaslin2', all_results) & grepl('^ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  fit_out_joint$feature <- gsub('.*s__', '', fit_out_joint$feature)
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_single
    fit_out_joint$qval <- fit_out_joint$qval_single
    if (grepl('Maaslin3ItAug', result)) {
      fit_out_joint$tool <- 'Maaslin3ItAug'
    } else {
      fit_out_joint$tool <- 'Maaslin3'
    }
  } else if (grepl('Maaslin2', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval
    fit_out_joint$qval <- fit_out_joint$qval
    fit_out_joint$tool <- 'Maaslin2'
    fit_out_joint$association <- 'abundance'
    fit_out_joint2 <- fit_out_joint
    fit_out_joint2$association <- 'prevalence'
    fit_out_joint <- rbind(fit_out_joint, fit_out_joint2)
  }
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool', 'association')]
  
  growing_df <- rbind(growing_df, fit_out_joint)
}

growing_df <- growing_df[growing_df$feature %in% keep_taxa,]
growing_df$association_name <- paste0(growing_df$feature, '-', growing_df$metadata_value, '-', growing_df$association)

for_scatter <- growing_df[growing_df$association %in% c('abundance'), c("association_name", "coef", "qval", "tool")]
split1 <- for_scatter[for_scatter$tool == 'Maaslin3',]
split2 <- for_scatter[for_scatter$tool == 'Maaslin3ItAug',]
split3 <- for_scatter[for_scatter$tool == 'Maaslin2',]
colnames(split1) <- c('association_name', 'Maaslin3', 'qval_Maaslin3', 'tool')
colnames(split2) <- c('association_name', 'Maaslin3ItAug', 'qval_Maaslin3ItAug', 'tool')
colnames(split3) <- c('association_name', 'Maaslin2', 'qval_Maaslin2', 'tool')
split1$tool <- NULL
split2$tool <- NULL
split3$tool <- NULL
joined_df <- full_join(split2, split1, by=c('association_name'))
joined_df <- full_join(joined_df, split3, by=c('association_name'))
joined_df$color <- paste0(ifelse(!is.na(joined_df$qval_Maaslin3ItAug) & joined_df$qval_Maaslin3ItAug < 0.1, 'v3ItAug,', ''),
                          ifelse(!is.na(joined_df$qval_Maaslin3) & joined_df$qval_Maaslin3 < 0.1, 'v3,', ''),
                          ifelse(!is.na(joined_df$qval_Maaslin2) & joined_df$qval_Maaslin2 < 0.1, 'v2', ''))
joined_df$association <- 'abundance'
joined_df_abun <- joined_df[!is.na(joined_df$color) & joined_df$color != '',]

for_scatter <- growing_df[growing_df$association %in% c('prevalence'), c("association_name", "coef", "qval", "tool")]
split1 <- for_scatter[for_scatter$tool == 'Maaslin3',]
split2 <- for_scatter[for_scatter$tool == 'Maaslin3ItAug',]
split3 <- for_scatter[for_scatter$tool == 'Maaslin2',]
colnames(split1) <- c('association_name', 'Maaslin3', 'qval_Maaslin3', 'tool')
colnames(split2) <- c('association_name', 'Maaslin3ItAug', 'qval_Maaslin3ItAug', 'tool')
colnames(split3) <- c('association_name', 'Maaslin2', 'qval_Maaslin2', 'tool')
split1$tool <- NULL
split2$tool <- NULL
split3$tool <- NULL
joined_df <- full_join(split2, split1, by=c('association_name'))
joined_df <- full_join(joined_df, split3, by=c('association_name'))
joined_df$color <- paste0(ifelse(!is.na(joined_df$qval_Maaslin3ItAug) & joined_df$qval_Maaslin3ItAug < 0.1, 'v3ItAug,', ''),
                          ifelse(!is.na(joined_df$qval_Maaslin3) & joined_df$qval_Maaslin3 < 0.1, 'v3,', ''),
                          ifelse(!is.na(joined_df$qval_Maaslin2) & joined_df$qval_Maaslin2 < 0.1, 'v2', ''))
joined_df$association <- 'prevalence'
joined_df_prev <- joined_df[!is.na(joined_df$color) & joined_df$color != '',]

joined_df <- rbind(joined_df_abun, joined_df_prev)
joined_df$association_name <- gsub('-|abundance|prevalence', ' ', joined_df$association_name)
joined_df <- joined_df[,c('association_name', 'Maaslin3ItAug', 'Maaslin3', 'Maaslin2', 'qval_Maaslin3ItAug', 'qval_Maaslin3', 'qval_Maaslin2', 'color', 'association')]

taxa_abun_vec <- rowMeans(taxa_table > 0.001)
names(taxa_abun_vec) <- gsub('.*s__', '', names(taxa_abun_vec))
taxa_abun_df <- data.frame('taxa' = names(taxa_abun_vec), 'prevalence' = taxa_abun_vec)
joined_df$taxa <- sapply(strsplit(joined_df$association_name, ' '), `[[`, 1)
joined_df <- left_join(joined_df, taxa_abun_df, by = 'taxa')
joined_df$taxa <- NULL

write.table(joined_df[rowSums(is.na(joined_df)) == 0,], file = 'tables/IBD.tsv', sep='\t', row.names = F)

##############
# HALLA prep #
##############

taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', sep = '\t', skip = 1)

colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
taxa_table <- taxa_table[taxa_table$clade_name == 'UNCLASSIFIED' | 
                           (grepl('\\|s__', taxa_table$clade_name) &
                              !grepl('\\|t__', taxa_table$clade_name)),]
rownames(taxa_table) <- taxa_table$clade_name
taxa_table$clade_name <- NULL
keep_taxa <- rownames(taxa_table)[rowMeans(taxa_table > 0.1) > 0.01]

# In-line number
mean(taxa_table[rowMeans(taxa_table > 0.1) > 0.01,] > 0)

fit_out_joint <- read.csv(paste0('results/ibd_associations_Maaslin3ItAug.tsv'), sep = '\t')
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[fit_out_joint$feature %in% make.names(keep_taxa),]

fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
fit_out_joint <- fit_out_joint[fit_out_joint$value %in% c('dysbiosis_CD'),]
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
fit_out_joint$pval <- fit_out_joint$pval_single
fit_out_joint$qval <- fit_out_joint$qval_single

fit_out_joint <- fit_out_joint[fit_out_joint$qval < 0.1,]

taxa_table_halla <- taxa_table[make.names(rownames(taxa_table)) %in% fit_out_joint$feature,]
taxa_table_halla <- taxa_table_halla[,order(colnames(taxa_table_halla))]
taxa_table_halla_colnames <- gsub('_P$', '', colnames(taxa_table_halla))
taxa_table_halla <- taxa_table_halla[,!duplicated(taxa_table_halla_colnames)]

mbx_table <- read.csv("data/intensities_hmp2.csv")
taxa_table_halla <- taxa_table_halla[,colnames(taxa_table_halla) %in% colnames(mbx_table)]
rownames(taxa_table_halla) <- gsub("_", " ", gsub(".*s__", "", rownames(taxa_table_halla)))

write.table(taxa_table_halla, 'data/taxa_table_halla.tsv', sep = '\t', row.names = T)







