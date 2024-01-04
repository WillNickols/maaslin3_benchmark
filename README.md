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