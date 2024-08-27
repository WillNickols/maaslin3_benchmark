# maaslin3_benchmark

## Installation

The necessary packages for running the benchmarking can be installed from the `environment.yml` file. After creating a Conda environment with these packages, run the following installations:

```
install.packages(c('dplyr', 'pbapply', 'lmerTest', 'parallel', 'lme4', 'plyr', 'optparse', 'logging', 'data.table', 'ggplot2', 'grid', 'pheatmap'))
install.packages(c("pkgmaker", "stringi", "doParallel", "SimSeq", "tidyr", "devtools", "TcGSA", "MCMCprecision")) # Come back to install devtools if necessary
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "SparseDOSSA2", "ALDEx2", "ANCOMBC", "TreeSummarizedExperiment", "Maaslin2"))
library("devtools")
install_github("biobakery/MaAsLin3")
```

## Evaluations

Each evaluation folder (`community_shift`, `deep_sequencing`, `general_evaluations`, `groups`, `omps`, and `unscaled`) is structured in the same way. Each evaluation has:
- A `data_generation` folder with scripts to generate data according to SparseDOSSA2 or ANCOM-BC's generator
- A `run_tools` folder with scripts to run each differential abundance tool on the generated datasets
- An `evaluate_tools` folder with scripts to compute accuracy metrics for each tool and produce plots
- A `.py` file with a workflow to generate the simulated data for all sets of parameters and run the differential abundance tools on all generated datasets
- A `.txt` file with the set of generators to use

The `library` folder contains general-purpose functions that apply over multiple evaluation types such as functions for data generation and accuracy evaluation.

### Running the evaluations

The following command runs the community shift evaluations.
```
python community_shift/evaluate_community_shift.py \
  --parameters community_shift/evaluate_community_shift.txt  \
  -o community_shift/
```

The following command runs the deep sequencing evaluations.
```
python deep_sequencing/evaluate_deep_sequencing.py \
  --parameters deep_sequencing/evaluate_deep_sequencing.txt  \
  -o deep_sequencing/
```

The following command runs the general evaluation with all generators.
```
python general_evaluations/evaluate_general.py \
  --parameters general_evaluations/evaluate_general.txt \
  -o general_evaluations/
```

The following command runs the groupwise difference evaluations.
```
python groups/evaluate_groups.py \
  --parameters groups/evaluate_groups.txt \
  -o groups/
```

The following command runs the ordered predictor evaluations.
```
python omps/evaluate_omps.py \
  --parameters omps/evaluate_omps.txt \
  -o omps/
```

The following command runs the spike-in abundance evaluations.
```
python unscaled/evaluate_unscaled.py \
  --parameters unscaled/evaluate_unscaled.txt  \
  -o unscaled/
```

## Randomization test

The `randomization_test` directory is structured similarly. It contains:
- A `run_tools` folder with scripts to run the randomized or non-randomized datasets
- An `evaluate_tools` folder with a script to compute accuracy metrics for each tool and produce plots
- A `.py` file with a workflow to generate the shuffled data and run the differential abundance tools on all generated datasets

The following command runs the randomization test evaluations.
```
python randomization_test/evaluate_randomization.py \
  -o randomization_test/
```

## Real data absolute abundance

The `real_data_absolute_abundance` directory contains three sub-directories, one for each dataset. Each contains:
- A `data` folder with the abundance data and metadata. These data were obtained from the supplementary information of each study or from the ENA nucleotide browser's display of per-sample metadata including read depth.
- A `run_scripts` folder with scripts to run the differential abundance tools on each dataset.
- A `join_results.R` script to combine the results and create plots.

## iHMP analysis

The `scripts` folder contains the script to perform the MetaPhlAn analysis of the iHMP data. The `results` folder contains the following:
- A `data` folder with the taxonomic profiles, metabolomic profiles, and patient metadata. Because of its size, the metabolomics file should be separately downloaded into this folder as `intensities_hmp2.csv`.
- A `run_scripts` folder with scripts to run each differential abundance tool
- A `diet_associations.py` script to run MaAsLin 3 for diet associations
- An `ibd_associations.py` script to run all differential abundance tools
- A `join_results.R` script to compile the taxonomic abundance (IBD and diet) results and create figures
- A `join_results_mbx.R` script to compile the metabolomics results and create figures

### Running the analysis

The following command performs the MetaPhlAn analysis.
```
python run_mpa.py -i data/hmp2_qc/ \
  -o maaslin3_benchmark/HMP2/outputs  \
  --input-extension fastq.gz --bowtie metaphlan4/
```

The following commands run the iHMP analysis.
```
python ibd_associations.py \
  -o maaslin3_benchmark/HMP2/analysis/ \
  --workingDirectory $(pwd)

python diet_associations.py \
  -omaaslin3_benchmark/HMP2/analysis_diet/ \
  --workingDirectory $(pwd)
```
