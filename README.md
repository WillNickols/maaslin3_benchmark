# maaslin3_benchmark

## Installation

The necessary packages for running the benchmarking can be installed from the `environment.yml` file. After creating a Conda environment with these packages, run the following installations:

```
install.packages(c('dplyr', 'pbapply', 'lmerTest', 'parallel', 'lme4', 'plyr', 'optparse', 'logging', 'data.table', 'ggplot2', 'grid', 'pheatmap'))
install.packages(c("pkgmaker", "stringi", "doParallel", "SimSeq", "tidyr", "devtools", "TcGSA", "MCMCprecision")) # Come back to install devtools if necessary
install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "SparseDOSSA2", "ALDEx2", "ANCOMBC", "TreeSummarizedExperiment", "Maaslin2"))
```

## Evaluations

Each evaluation folder (`general_evaluations`, `groups`, `omps`, and `unscaled`) is structured in the same way. Each evaluation has:
- A `data_generation` folder with scripts to generate data according to SparseDOSSA2, ANCOM-BC's generator, or SimSeq as applicable
- A `run_tools` folder with scripts to run each differential abundance tool on the generated datasets
- An `evaluate_tools` folder with scripts to compute accuracy metrics for each tool. (The `evaluate_tools` folder in the `general_evaluations` folder also contains the code to generate the thesis figures.)
- A `.py` file with a workflow to generate the simulated data for all sets of parameters and run the differential abundance tools on all generated datasets
- A `.txt` file with the set of generators to use

The `library` folder also contains general-purpose functions that apply over multiple evaluation types.

### Running the evaluations

The following command runs the general evaluation with all generators.
```
python general_evaluations/evaluate_general.py \
  --parameters general_evaluations/evaluate_general.txt \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/general_evaluations/ \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 12 --time 1200 --mem 30000 \
  --local-jobs 12
```

The following command runs the OMP and GOMP evaluations.
```
python omps/evaluate_omps.py \
  --parameters omps/evaluate_omps.txt \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/omps/ \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 20 --time 240 --mem 40000 \
  --local-jobs 12
```

The following command runs the group (categorical variable) evaluations.
```
python groups/evaluate_groups.py \
  --parameters groups/evaluate_groups.txt \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/groups/ \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 8 --time 240 --mem 10000 \
  --local-jobs 12
```

The following command runs the spike-in abundance evaluations.
```
python unscaled/evaluate_unscaled.py \
  --parameters unscaled/evaluate_unscaled.txt  \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/unscaled/  \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 12 --time 240 --mem 30000  \
  --local-jobs 12
```

## iHMP analysis

The `scripts` folder contains the script to perform the MetaPhlAn analysis of the iHMP data and the command to perform the HAllA analysis. The `results` folder contains the following:
- A `data` folder with the taxonomic profiles, metabolomic profiles, and patient metadata. Because of its size, the metabolomics file should be separately downloaded into this folder as `intensities_hmp2.csv`.
- A `run_scripts` folder with scripts to run each differential abundance tool
- A `diet_associations.py` script to run MaAsLin 3 for diet associations
- An `ibd_associations.py` script to run all differential abundance tools
- A `join_results.R` script to compile the taxonomic abundance (IBD and diet) results and create figures
- A `join_results_mbx.R` script to compile the metabolomics results and create figures

### Running the analysis

The following command performs the MetaPhlAn analysis.
```
python run_mpa.py -i /n/hutlab12_nobackup/data/hmp2_qc/ \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/HMP2/outputs  \
  --grid-partition 'shared' --grid-jobs 200 --cores 15 --time 240 --mem 30000 \
  --input-extension fastq.gz --bowtie /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/metaphlan4/ \
  --local-jobs 12 --grid-scratch /n/holyscratch01/huttenhower_lab/wnickols/HMP2/
```

The following commands run the iHMP analysis.
```
python ibd_associations.py \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/HMP2/analysis/ \
  --workingDirectory /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/ \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 12 --time 240 --mem 80000   --local-jobs 12

python diet_associations.py \
  -o /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/maaslin3_benchmark/HMP2/analysis_diet/ \
  --workingDirectory /n/holystore01/LABS/huttenhower_lab/Users/wnickols/maaslin3/ \
  --grid-partition 'sapphire' --grid-jobs 96 --cores 12 --time 240 --mem 80000   --local-jobs 12
```

[//]: # ( HALLA installed by creating conda environment for rpy2, pip installing all the other packages, editing out the sklearn requirement of requirements.txt and using setup.py )



