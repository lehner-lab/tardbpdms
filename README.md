# Overview

Welcome to the GitHub repository for the following publication: "The mutational landscape of a prion-like domain" (Bolognesi B & Faure AJ et al., 2019)

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Required Software

To run the tardbpdms pipeline you will need the following software and associated packages:

* **[R](https://www.r-project.org/) >=v3.5.2** (Biostrings, caTools, corpcor, cowplot, data.table, gdata, ggplot2, GGally, hexbin, lemon, optparse, parallel, pdist, plyr, ppcor, raster, reshape2, Rpdb, RColorBrewer)

The following packages are optional:

* **[DiMSum](https://github.com/lehner-lab/DiMSum)** (pipeline for pre-processing deep mutational scanning data i.e. FASTQ to counts)
* **[DMS2structure](https://github.com/lehner-lab/DMS2structure)** (scripts used for epistasis and structure analysis of deep mutational scanning data in Schmiedel & Lehner, bioRxiv 2018)

# Installation and loading

Open R and enter:

```
# Install
if(!require(devtools)) install.packages("devtools")
devtools::install_github("lehner-lab/tardbpdms")

# Load
library(tardbpdms)

# Help
?tardbpdms
```

# Required Data

All pre-processed data and required miscellaneous files should be downloaded from [here](https://www.dropbox.com/s/nbuvtfler395rcd/tardbpdms_misc.zip?dl=0) to your base (project) directory i.e. where output files should be written, and unzipped.

# Pipeline

The top-level function **tardbpdms()** is the recommended entry point to the pipeline and reproduces the figures and results from the computational analyses described in the following publication: "The mutational landscape of a prion-like domain" (Bolognesi B & Faure AJ et al., 2019). See section on "Required Data" above for instructions on how to obtain all required data and miscellaneous file before running the pipeline.

## Stage 1: DiMSum counts to fitness

**tardbpdms_dimsumcounts_to_fitness**

## Stage 2: Quality control plots

**tardbpdms_quality_control**

## Stage 3: Combine toxicity estimates from 290 and 332 experiments

**tardbpdms_combine_toxicity**

## Stage 4: Calculate single and double mutant effects from AA PCA

**tardbpdms_aa_properties_mutant_effects**

## Stage 5: Calculate single and double mutant effects from aggregation tool predictions

**tardbpdms_agg_tools_mutant_effects**

## Stage 6: Single mutant heatmaps

**tardbpdms_single_mutant_heatmaps**

## Stage 7: Human disease mutations

**tardbpdms_human_disease_mutations**

## Stage 8: Dot plots showing explained variance of models to predict variant toxicity

**tardbpdms_toxicity_model_summary**

## Stage 9: Violin plots showing toxicity vs #introduced AAs

**tardbpdms_num_introduced_aa_violins**

## Stage 10: Hydrophobicity of WT and toxicity hotspot line plot

**tardbpdms_wt_hydrophobicity**

## Stage 11: Epistasis analysis

**tardbpdms_epistasis_analysis**

## Stage 12: Secondary structure predictions

**tardbpdms_secondary_structure_predictions**

## Stage 13: Guenther structure propensities

**tardbpdms_guenther_structure_propensities**

## Stage 14: PWI heatmaps

**tardbpdms_PWI_heatmaps**

