library(UpSetR)
library(reshape2)
library(dplyr)
library(ggplot2)
library(stringr)

workingDirectory <- "~/Documents/GitHub/maaslin3_benchmark"
figures_folder <- paste0(workingDirectory, '/Figures/thesis_figures/')

Maaslin3_path <- paste0(gsub("/$", "", workingDirectory), "/Maaslin3/R/")
for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  source(file.path(Maaslin3_path, R_file))
}

# Read in data
mbx_table <- read.csv("data/intensities_hmp2.csv")
annotations <- read.csv("data/annotations_hmp2.csv")

if (all(abs(annotations$prev - rowMeans(!is.na(mbx_table))) < 0.001)) {
  mbx_table <- mbx_table[annotations$prim_feature == 'primary' & annotations$Metabolite != '',]
  rownames(mbx_table) <- paste0(annotations$Metabolite[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                annotations$HMDB.ID[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                annotations$Method[annotations$prim_feature == 'primary' & annotations$Metabolite != ''])
}

keep_taxa <- rownames(mbx_table)[!grepl('redundant', rownames(mbx_table))]

all_results <- list.files('results/')

growing_df <- data.frame()
for (result in all_results[grepl('^mbx_ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]

  if (grepl('Maaslin2', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
    fit_out_joint$tool <- 'Maaslin2'
  }
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_joint
    fit_out_joint$qval <- fit_out_joint$qval_joint
    fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
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

growing_df$feature <- make.names(growing_df$feature)
keep_taxa <- make.names(keep_taxa)
growing_df <- growing_df[growing_df$qval < 0.1 & growing_df$feature %in% keep_taxa,]
growing_df$association <- paste0(growing_df$feature, '-', growing_df$metadata_value)

listInput <- list(`Maaslin 2` = growing_df$association[growing_df$tool == 'Maaslin2'], 
     `Maaslin 3` = growing_df$association[growing_df$tool == 'Maaslin3'], 
     `Maaslin 3 Iterative Augmented` = growing_df$association[growing_df$tool == 'Maaslin3ItAug'], 
     `ANCOM-BC2` = growing_df$association[growing_df$tool == 'ANCOMBC'],
     ALDEx2 = growing_df$association[growing_df$tool == 'ALDEx2'])

plot_out <- upset(fromList(listInput), order.by = "freq", text.scale = 2.5)

# In-text number
length(unique(listInput$`Maaslin 2`) %in% unique(listInput$`Maaslin 3`)) / length(unique(c(listInput$`Maaslin 2`, listInput$`Maaslin 3`)))

png(file=paste0(figures_folder, 'fig_16a.png'), width = 7, height = 6, res = 1000, units = 'in')
plot_out
dev.off()

# Scatter plot
growing_df <- data.frame()
for (result in all_results[grepl('Maaslin3|Maaslin2', all_results) & grepl('^mbx_ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
  
  if (grepl('Maaslin2', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
    fit_out_joint$tool <- 'Maaslin2'
  }
  
  if (grepl('Maaslin3', result)) {
    fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
    fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
    fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
    fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
    fit_out_joint$pval <- fit_out_joint$pval_joint
    fit_out_joint$qval <- fit_out_joint$qval_joint
    fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
    if (grepl('Maaslin3ItAug', result)) {
      fit_out_joint$tool <- 'Maaslin3ItAug'
    } else {
      fit_out_joint$tool <- 'Maaslin3'
    }
    fit_out_joint <- fit_out_joint[fit_out_joint$association == 'abundance',]
  }
  
  fit_out_joint <- fit_out_joint[,c('feature', 'metadata_value', 'coef', 'pval', 'qval', 'tool')]
  
  growing_df <- rbind(growing_df, fit_out_joint)
}

growing_df$feature <- make.names(growing_df$feature)
growing_df <- growing_df[growing_df$feature %in% keep_taxa,]
growing_df$association_name <- paste0(growing_df$feature, '-', growing_df$metadata_value)

for_scatter <- growing_df[, c("association_name", "coef", "qval", "tool")]
split1 <- for_scatter[for_scatter$tool == 'Maaslin3',]
split2 <- for_scatter[for_scatter$tool == 'Maaslin2',]
colnames(split1) <- c('association_name', 'Maaslin3', 'qval_Maaslin3', 'tool')
colnames(split2) <- c('association_name', 'Maaslin2', 'qval_Maaslin2', 'tool')
split1$tool <- NULL
split2$tool <- NULL
joined_df <- left_join(split2, split1, by=c('association_name'))
joined_df$color <- ifelse(joined_df$qval_Maaslin2 < 0.1,
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Both Significant', 
                                 'Maaslin 2 Only'),
                          ifelse(joined_df$qval_Maaslin3 < 0.1, 'Maaslin 3 Only', 'FF'))
joined_df <- joined_df[!is.na(joined_df$color) & joined_df$color != 'FF',]

plot_out <- ggplot(joined_df, aes(x = Maaslin2, 
                      y = Maaslin3, color = color)) + 
  geom_point(size=3, alpha=0.5) + 
  theme_bw() + 
  theme(text = element_text(size = 25)) + 
  geom_abline(slope=1, intercept=0, linetype = 'dashed') + 
  labs(color = '', x = 'Maaslin 2 Coefficient',
       y = 'Maaslin 3 Coefficient') + 
  scale_color_brewer(palette="Set1")
ggsave(paste0(figures_folder, 'fig_16c.png'),
       plot = plot_out, width = 8, height = 5)

# Volcano plot
fit_out_joint <- read.csv(paste0('results/mbx_ibd_associations_Maaslin3.tsv'), sep = '\t')
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[fit_out_joint$feature %in% keep_taxa,]

fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
fit_out_joint$pval <- fit_out_joint$pval_single
fit_out_joint$qval <- fit_out_joint$qval_single

significance_threshold <- 0.1
fit_out_joint$color <- ifelse(fit_out_joint$qval < significance_threshold, 
                              ifelse(fit_out_joint$association == 'prevalence', 'Significant prevalence', 'Significant abundance'),
                              ifelse(fit_out_joint$association == 'prevalence', 'Insignificant prevalence', 'Insignificant abundance'))
qval_threshold <- 20
fit_out_joint$qval <- ifelse(fit_out_joint$qval < 10^-qval_threshold, 0, fit_out_joint$qval)
fit_out_joint$coef <- ifelse(fit_out_joint$coef > 25, Inf, fit_out_joint$coef)
volcano_plot <- ggplot(fit_out_joint[sample(1:nrow(fit_out_joint)),], aes(x = coef, y = -log10(qval))) +
  geom_point(aes(color = as.factor(color)), size = 4, alpha=0.5) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
  labs(x = "Effect Size", y = "-log10(q-value)") +
  theme_bw() + 
  theme(text = element_text(size = 35),
        legend.position = 'right',
        legend.direction = 'vertical') +
  labs(color = '') + 
  ylim(c(0, qval_threshold)) +
  xlim(c(-25, 25)) +
  scale_color_manual(values=c("#ff9496", "#9bd0e6", "#f71c1a", "#007ab9"))

# In-text number
table(fit_out_joint$association, fit_out_joint$coef > 0) / sum(table(fit_out_joint$association))

ggsave(paste0(figures_folder, 'fig_16d.png'),
       plot = volcano_plot, width = 18, height = 7.2)

##########################
# Abundance associations #
##########################

fit_out_joint <- read.csv(paste0('results/mbx_ibd_associations_Maaslin3.tsv'), sep = '\t')
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[fit_out_joint$feature %in% keep_taxa,]

fit_out_joint <- fit_out_joint[fit_out_joint$metadata %in% c('diagnosis', 'dysbiosis_state'),]
fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]
fit_out_joint <- fit_out_joint[is.na(fit_out_joint$error),]
fit_out_joint$metadata_value <- paste0(fit_out_joint$metadata, '_', fit_out_joint$value)
fit_out_joint$pval <- fit_out_joint$pval_single
fit_out_joint$qval <- fit_out_joint$qval_single
fit_out_joint$hmdb <- gsub("\\..*|_.*", "", gsub(".*_(HMDB)", "\\1", fit_out_joint$feature))
fit_out_joint <- fit_out_joint[grepl('^HMDB', fit_out_joint$hmdb),]

HMDB_chemical_taxonomy <- read.csv('data/HMDB_chemical_taxonomy.csv')
fit_out_joint <- left_join(fit_out_joint, HMDB_chemical_taxonomy, by=c('hmdb'='HMDB.ID'))
fit_out_joint <- fit_out_joint[fit_out_joint$value == 'dysbiosis_CD',]
#fit_out_joint <- fit_out_joint[fit_out_joint$qval_single < 0.1,]

fit_out_joint$org_feature <- fit_out_joint$feature
fit_out_joint$feature <- gsub('_.*', '', fit_out_joint$feature)
fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
fit_out_joint$feature <- gsub("(\\d)\\.(\\d)", "\\1,\\2", fit_out_joint$feature)
fit_out_joint$feature <- gsub("(\\d)\\.", "\\1-", fit_out_joint$feature)
fit_out_joint$feature <- gsub("\\.(\\d)", "-\\1", fit_out_joint$feature)
fit_out_joint$feature <- gsub("\\.", " ", fit_out_joint$feature)
fit_out_joint$feature <- gsub("(C[0-9]+),([0-9]+)", "\\1:\\2", fit_out_joint$feature)
fit_out_joint$feature <- gsub("(C[0-9,:]+)-", "\\1 ", fit_out_joint$feature)
fit_out_joint$feature <- gsub("([0-9])- ([A-Z]) ", "\\1'-\\2-", fit_out_joint$feature)

for (sub_class in unique(fit_out_joint$Sub_Class)) {
  fit_out_joint$subclass_signif[fit_out_joint$Sub_Class == sub_class] <- 
    wilcox.test(fit_out_joint[fit_out_joint$Sub_Class == sub_class,]$coef)$p.value
}

fit_out_joint$subclass_signif <- mapvalues(fit_out_joint$subclass_signif, 
                                           unique(fit_out_joint$subclass_signif), 
                                           p.adjust(unique(fit_out_joint$subclass_signif)))

fit_out_joint <- fit_out_joint[fit_out_joint$association == 'abundance' & 
                fit_out_joint$Sub_Class != '' &
                fit_out_joint$Sub_Class %in% names(which(table(fit_out_joint$Sub_Class) > 5)),]
fit_out_joint$Sub_Class <- factor(fit_out_joint$Sub_Class, levels = sort(unique(fit_out_joint$Sub_Class), decreasing = T))

fit_out_joint$signif_star <- case_when(fit_out_joint$subclass_signif < 0.001 ~ "***",
                                       fit_out_joint$subclass_signif < 0.01 ~ "**",
                                       fit_out_joint$subclass_signif < 0.1 ~ "*",
                                       fit_out_joint$subclass_signif <= 1 ~ "")
fit_out_joint <- fit_out_joint[order(fit_out_joint$Sub_Class, decreasing = F),]
labs1 <- levels(fit_out_joint$Sub_Class)
labs2 <- fit_out_joint$signif_star[!duplicated(fit_out_joint$Sub_Class)]

plot_out <- ggplot(fit_out_joint, 
       aes(x = as.numeric(Sub_Class), y = coef)) + 
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

ggsave(paste0(figures_folder, 'fig_17a.png'),
       plot = plot_out, width = 14, height = 6)
  
top_hits <- fit_out_joint[fit_out_joint$association == 'abundance' & fit_out_joint$feature != 'UNCLASSIFIED',]
top_hits <- rbind(top_hits[top_hits$coef > 0,][order(top_hits[top_hits$coef > 0,]$qval_single)[1:8],],
                  top_hits[top_hits$coef < 0,][order(top_hits[top_hits$coef < 0,]$qval_single)[1:8],])
top_hits <- rbind(top_hits[order(-top_hits$coef)[1:8],],
                  top_hits[order(top_hits$coef)[8:1],])
top_hits$feature <- gsub('_', ' ', top_hits$feature)
top_hits$feature <- paste0(str_to_title(top_hits$feature))
top_hits$feature <- paste0(top_hits$feature, ' (', top_hits$Sub_Class, ')')
top_hits$feature <- factor(top_hits$feature, levels = unique(top_hits$feature))

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

ggsave(paste0(figures_folder, 'fig_17b.png'),
       plot = plot_out, width = 14, height = 6)

#################
# Results table #
#################

growing_df <- data.frame()
for (result in all_results[grepl('Maaslin3|Maaslin2', all_results) & grepl('^mbx_ibd_', all_results)]) {
  fit_out_joint <- read.csv(paste0('results/', result), sep = '\t')
  
  fit_out_joint <- fit_out_joint[!is.na(fit_out_joint$feature),]

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
  fit_out_joint$feature <- gsub("^X([0-9]+)", "\\1", fit_out_joint$feature)
  
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

taxa_abun_vec <- rowMeans(mbx_table)
names(taxa_abun_vec) <- make.names(names(taxa_abun_vec))
taxa_abun_df <- data.frame('taxa' = names(taxa_abun_vec), 'abundance' = taxa_abun_vec)
joined_df$taxa <- sapply(strsplit(joined_df$association_name, ' '), `[[`, 1)
joined_df <- left_join(joined_df, taxa_abun_df, by = 'taxa')
joined_df$taxa <- NULL

write.table(joined_df[rowSums(is.na(joined_df)) == 0,], file = 'tables/mbx_IBD.tsv', sep='\t', row.names = F)

#########
# HALLA #
#########

mbx_table <- read.csv("data/intensities_hmp2.csv")
annotations <- read.csv("data/annotations_hmp2.csv")

if (all(abs(annotations$prev - rowMeans(!is.na(mbx_table))) < 0.001)) {
  mbx_table <- mbx_table[annotations$prim_feature == 'primary' & annotations$Metabolite != '',]
  rownames(mbx_table) <- paste0(annotations$Metabolite[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                annotations$HMDB.ID[annotations$prim_feature == 'primary' & annotations$Metabolite != ''], '_',
                                annotations$Method[annotations$prim_feature == 'primary' & annotations$Metabolite != ''])
}
keep_taxa <- rownames(mbx_table)[!grepl('redundant', rownames(mbx_table))]

fit_out_joint <- read.csv(paste0('results/mbx_ibd_associations_Maaslin3.tsv'), sep = '\t')
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

mbx_table_halla <- mbx_table[make.names(rownames(mbx_table)) %in% fit_out_joint$feature,]
mbx_table_halla[is.na(mbx_table_halla)] <- 0

taxa_table <- read.csv('data/metaphlan4_taxonomic_profiles.tsv', sep = '\t', skip = 1)
colnames(taxa_table) <- gsub('_taxonomic$', '', colnames(taxa_table))
mbx_table_halla <- mbx_table_halla[,colnames(mbx_table_halla) %in% colnames(taxa_table)]

mbx_table_halla_rownames <- gsub('_.*', '', rownames(mbx_table_halla))
mbx_table_halla_rownames <- gsub("^X([0-9]+)", "\\1", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("(\\d)\\.(\\d)", "\\1,\\2", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("(\\d)\\.", "\\1-", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("\\.(\\d)", "-\\1", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("\\.", " ", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("(C[0-9]+),([0-9]+)", "\\1:\\2", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("(C[0-9,:]+)-", "\\1 ", mbx_table_halla_rownames)
mbx_table_halla_rownames <- gsub("([0-9])- ([A-Z]) ", "\\1'-\\2-", mbx_table_halla_rownames)
mbx_table_halla_rownames <- paste(toupper(substr(mbx_table_halla_rownames, 1, 1)), substr(mbx_table_halla_rownames, 2, nchar(mbx_table_halla_rownames)), sep="")

mbx_table_halla <- mbx_table_halla[!duplicated(mbx_table_halla_rownames),]
rownames(mbx_table_halla) <- mbx_table_halla_rownames[!duplicated(mbx_table_halla_rownames)]

write.table(mbx_table_halla, 'data/mbx_table_halla.tsv', sep = '\t', row.names = T)




