remove(list = ls())

dataset_names <- c(ArcticFireSoils="Soil - Fires",
                   ArcticFreshwaters="Freshwater - Arctic",
                   ArcticTransects="Soil - Arctic",
                   art_scher="Human - RA",
                   asd_son= "Human - ASD",
                   BISCUIT= "Human - CD (1)",
                   Blueberry= "Soil - Blueberry",
                   cdi_schubert="Human - C. diff (1)",
                   cdi_vincent="Human - C. diff (2)",
                   Chemerin="Mouse Facilities",
                   crc_baxter="Human - CC (1)",
                   crc_zeller="Human - CC (2)",
                   edd_singh="Human - Inf.",
                   Exercise="Mouse - Exercised",
                   glass_plastic_oberbeckmann="Marine - Plastic (4)",
                   GWMC_ASIA_NA="WWSR - Continents",
                   GWMC_HOT_COLD="WWSR - Temp.",
                   hiv_dinh="Human - HIV (1)",
                   hiv_lozupone="Human - HIV (2)",
                   hiv_noguerajulian="Human - HIV (3)",
                   ibd_gevers="Human - CD (2)",
                   ibd_papa="Human - IBD",
                   Ji_WTP_DS="Freshwater - Treat.",
                   MALL="Human - ALL",
                   ob_goodrich="Human - OB (1)",
                   ob_ross="Human - OB (2)",
                   ob_turnbaugh="Human - OB (3)",
                   ob_zhu="Human - OB (4)",
                   Office="Built - Office",
                   par_scheperjans="Human - Par.",
                   sed_plastic_hoellein="Marine - Plastic (2)",
                   sed_plastic_rosato="Marine - Plastic (5)",
                   seston_plastic_mccormick="River - Plastic",
                   sw_plastic_frere="Marine - Plastic (1)",
                   sw_sed_detender="Marine - Sediment",
                   t1d_alkanani="Human - T1D (1)",
                   t1d_mejialeon="Human - T1D (2)",
                   wood_plastic_kesy="Marine - Plastic (3)")

included_datasets <- c("Built - Office",
                       "Freshwater - Arctic",
                       "Freshwater - Treat.",
                       "Human - C. diff (1)",
                       "Human - HIV (3)",
                       "Human - OB (1)",
                       "Marine - Sediment",
                       "Soil - Blueberry")

package_vec = c("reshape2", "ggplot2", "optparse", 
                "parallel", "stringi", "parallel", "plyr", "tidyr", "scales", "pbapply")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))
source('library/run_evaluation_helpers.R')
library(ggbeeswarm)

files_in <- list.files('randomization_test/associations/', full.names = T)

signif_threshold <- 0.05

process_file <- function(file) {
  in_data <- read.csv(file, sep = '\t')
  
  prop_below_threshold <- in_data %>%
    dplyr::group_by(dataset, associations) %>%
    dplyr::summarize(
      qval = mean(qval < signif_threshold, na.rm = TRUE),
      qval_joint = mean(qval_joint < signif_threshold, na.rm = TRUE)
    )
  
  if (grepl("non_null.tsv", file)) {
    prop_below_threshold$iter <- "non_null"
  } else {
    prop_below_threshold$iter <- gsub(".*_|\\.tsv", "", file)
  }
  
  return(prop_below_threshold)
}

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
clusterEvalQ(cl, library(dplyr))
clusterExport(cl, c("process_file", "signif_threshold"))

result_list <- pblapply(files_in, process_file, cl = cl)
growing_df <- do.call(rbind, result_list)

stopCluster(cl)

results_df <- reshape2::melt(growing_df, id.vars = c('dataset', 'iter', "associations"))
colnames(results_df) <- c("dataset", "iter", "associations", "qval_type", "value")
results_df <- reshape2::melt(results_df, id.vars = c('dataset', 'iter', 'qval_type', 'associations'))
results_df <- unique(results_df)
results_df <- results_df %>%
  mutate(qval_type = case_when(
    qval_type == 'qval_joint' & associations == 'abundance' ~ 'Joint Q-value',
    qval_type == 'qval' & associations == 'abundance' ~ 'Abundance Q-value',
    qval_type == 'qval' & associations == 'prevalence' ~ 'Prevalence Q-value'
  )
  )
results_df <- results_df %>%
  mutate(qval_type = factor(qval_type, levels = c('Joint Q-value', 'Abundance Q-value', 'Prevalence Q-value')))
results_df <- results_df[!is.na(results_df$qval_type),]

results_df$dataset <- dataset_names[results_df$dataset]
results_df$value <- results_df$value * 100
results_df$Data <- ifelse(results_df$iter == 'non_null', "Original", "Randomized")
results_df <- results_df[!is.na(results_df$dataset),] # %in% included_datasets

# In-text numbers
results_df %>%
    dplyr::filter(Data == 'Randomized') %>%
    dplyr::group_by(qval_type) %>%
    dplyr::summarise(max(value))

# print(results_df %>%
#     dplyr::filter(Data == 'Randomized') %>%
#     dplyr::group_by(qval_type, dataset) %>%
#     dplyr::summarise(mean(value == 0)), n = 1000)

# print(results_df %>%
#           dplyr::filter(Data == 'Randomized') %>%
#           dplyr::group_by(qval_type) %>%
#           dplyr::summarise(mean(value > 0)), n = 1000)

print(results_df %>%
          dplyr::filter(Data == 'Randomized') %>%
          dplyr::group_by(qval_type) %>%
          dplyr::summarise(sum(value > 0)), n = 1000)

print(results_df %>%
          dplyr::filter(Data == 'Original') %>%
          dplyr::group_by(qval_type, dataset) %>%
          dplyr::summarise(mean(value)), n = 1000) %>%
    dplyr::summarise(mean(`mean(value)`))

# mean((results_df %>%
#           dplyr::filter(Data == 'Original') %>%
#           dplyr::group_by(qval_type, dataset) %>%
#           dplyr::summarise(mean_value = mean(value)) %>%
#           tidyr::pivot_wider(names_from = qval_type, values_from = mean_value) %>%
#           dplyr::mutate(difference = `Prevalence Q-value` - `Abundance Q-value`))$difference)

plot_out <- ggplot(results_df, aes(x = dataset, y = value, color = Data)) + 
  geom_beeswarm(cex = 0.1, size = 2.5) + 
  facet_grid(. ~ qval_type) +
  labs(x = "Dataset",
    y = "Percent significant") +
  scale_y_continuous(breaks = c(0, 1, 2, 5, 10, 30, 100), 
                     trans = scales::pseudo_log_trans(base = 10)) +
  theme_bw() + 
  coord_flip() + 
  scale_color_manual(values = c('Original' = 'navy', 'Randomized' = 'firebrick')) + 
  theme(text=element_text(size=21),
        legend.position = 'right',
        strip.background = element_rect(fill = "gray95"))

figures_folder <- paste0('Figures/paper_figures/')
ggsave(paste0(figures_folder, 'randomization_test.png'),
       plot = plot_out, width = 12, height = 10)


