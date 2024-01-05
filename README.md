# maaslin3_benchmark

## General evaluations

### Data generation

One script for each of the following given parameters. These should produce a folder for each iteration of each set of parameters with `abundances.tsv`, `metadata.tsv`, and `truth.tsv`.
- SparseDOSSA2 on the expanded Maaslin 2 set of parameters
- ANCOM-BC's methods
- SimSeq on a few of Jacob's datasets

One workflow script for each of those to run all parameters.

### Per-tool outputs

One script for each of the following given the parameters. These should produce an `associations.tsv` file for each iteration of each set of parameters with the metadata, the taxon, the effect size, and the significance.
- Maaslin 2
- Maaslin 3
- Maaslin 3 (iterative, augmented)
- ANCOM-BC
- ALDEx2

One workflow script for each of these to run all parameter combinations.

### Results and visualizations

One script that produces a `results/` folder of everything.
- Metrics: sensitivity, specificity, FDR, F1 score, MCC, AUC, pAUC, conservative Area, liberal Area, total area, time
- What's being missed
- Effect sizes

## Maaslin 3 new inference evaluations

## Running the workflows

Installation: Conda environment, Maaslin3 packages, other R packages, Bioconductor
install.packages(c('dplyr', 'pbapply', 'lmerTest', 'parallel', 'lme4', 'plyr', 'optparse', 'logging', 'data.table', 'ggplot2', 'grid', 'pheatmap'))
install.packages(c("pkgmaker", "stringi", "doParallel", "SimSeq", "tidyr")) # Come back to install devtools if necessary
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "SparseDOSSA2", "ALDEx2", "ANCOMBC", "TreeSummarizedExperiment", "Maaslin2"))


Running:
```
python3 general_evaluations/data_generation/ANCOM_BC_generator_workflow.py --parameters general_evaluations/data_generation/ANCOM_BC_generator_tmp.txt --working-directory /Users/williamnickols/Documents/GitHub/maaslin3_benchmark --cores 4

python3 general_evaluations/data_generation/SD2_workflow.py --parameters general_evaluations/data_generation/SD2_tmp.txt --working-directory /Users/williamnickols/Documents/GitHub/maaslin3_benchmark --cores 4

python3 general_evaluations/data_generation/SimSeq_workflow.py --parameters general_evaluations/data_generation/SimSeq_tmp.txt --working-directory /Users/williamnickols/Documents/GitHub/maaslin3_benchmark --cores 4

python3 general_evaluations/run_tools/workflow.py --generators ANCOM_BC_generator,SD2,SimSeq --working-directory /Users/williamnickols/Documents/GitHub/maaslin3_benchmark --cores 4 --tmp
```

```
python general_evaluations/evaluate_general.py \
  --parameters general_evaluations/evaluate_general.txt \
  -o /n/hutlab12_nobackup/users/wnickols/maaslin3/maaslin3_benchmark \
  --grid-scratch /n/holyscratch01/huttenhower_lab/wnickols/maaslin3_benchmark/ \
  --grid-partition 'shared' --grid-jobs 96 --cores 8 --time 180 --mem 2000 \
  --local-jobs 12
```








