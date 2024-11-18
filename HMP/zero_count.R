library(curatedMetagenomicData)

input_data <- curatedMetagenomicData("HMP_2012.relative_abundance", dryrun = FALSE, rownames = "short")
mean((assay(input_data$`2021-03-31.HMP_2012.relative_abundance`) == 0))