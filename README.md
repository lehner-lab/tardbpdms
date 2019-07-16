# Overview

Welcome to the GitHub repository for the following publication: [The mutational landscape of a prion-like domain (Bolognesi B & Faure AJ et al., 2019)](https://www.biorxiv.org/content/10.1101/592121v1)

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

Variant counts, pre-processed data and required miscellaneous files should be downloaded from [here](https://www.dropbox.com/s/bnyjlsylak843qd/misc.zip?dl=0) to your project directory (see 'base_dir' argument) i.e. where output files should be written, and unzipped.

# Running

There are a number of options available for running the tardbpdms pipeline depending on user requirements.

* ## Basic (default)

Default pipeline functionality uses variant counts (see 'Required Data') to reproduce all figures in the publication. Neither **[DiMSum](https://github.com/lehner-lab/DiMSum)** nor **[DMS2structure](https://github.com/lehner-lab/DMS2structure)** packages are required for this default functionality.

* ## Raw read processing

Raw read processing is not handled by the tardbpdms pipeline. FastQ files ([GSE128165](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128165)) from paired-end sequencing of replicate deep mutational scanning (DMS) libraries before ('input') and after selection ('output') were processed using **[DiMSum](https://github.com/lehner-lab/DiMSum)** (manuscript in prep.), an R package that wraps common biological sequence processing tools.

DiMSum command-line arguments and Experimental design files required to obtain variant counts from FastQ files are available [here](https://www.dropbox.com/sh/dg609u5zalkozpn/AAAnJQKkR_cP5IaaOXOGGtApa?dl=0).

* ## Estimating fitness (1-toxicity) from variant counts

Pipeline stage 1 ('tardbpdms_dimsumcounts_to_fitness') estimates toxicity and error of single and double AA mutants from variant counts for each library separately. This stage is computationally intensive (~2hours on 10 cores) and is therefore not run by default. **Note:** When running the pipeline for the first time or to force re-execution of this stage set 'rerun_fitness = T'.

* ## Epistasis analysis

Pipeline stage 11 ('tardbpdms_epistasis_analysis') performs epistasis calculations separately for each DMS library. This stage is computationally intensive (~1hour on 10 cores) and is therefore not run by default. **Note:** 'Required Data' (see above) already includes precomputed results of the epistasis analysis necessary to reproduce the corresponding figures in the publication. However, to force re-execution of this stage set 'rerun_epistasis = T'. Additionally, the correct path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

* ## Structure analyses

Pipeline stages 12 and 13 ('tardbpdms_secondary_structure_predictions', 'tardbpdms_guenther_structure_propensities') perform secondary structure predictions and structure propensity calculations for PDB-structure derived
contact matrices respectively. Secondary structure predictions and propensity calculations are computationally intensive and are therefore not re-run by default. **Note:** 'Required Data' (see above) already includes precomputed results of the structure analyses necessary to reproduce the corresponding figures in the publication. To force re-execution set 'rerun_structure = T'. Additionally, the correct path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

# Pipeline

The top-level function **tardbpdms()** is the recommended entry point to the pipeline and reproduces the figures and results from the computational analyses described in the following publication: "The mutational landscape of a prion-like domain" (Bolognesi B & Faure AJ et al., 2019). See section on "Required Data" above for instructions on how to obtain all required data and miscellaneous files before running the pipeline.

## Stage 1: DiMSum counts to fitness

This stage ('tardbpdms_dimsumcounts_to_fitness') estimates toxicity and error of single and double AA mutants from variant counts for each library separately. This stage is computationally intensive (~2hours on 10 cores) and is therefore not run by default. When running the pipeline for the first time or to force re-execution of this stage set 'rerun_fitness = T'.

## Stage 2: Quality control plots

This stage ('tardbpdms_quality_control') produces quality control plots of toxicity estimates before and after inter-replicate normalisation.

## Stage 3: Combine toxicity estimates from 290 and 332 experiments

This stage ('tardbpdms_combine_toxicity') performs inter-library normalisation, toxicity distribution plots, growth rate comparison plots and position-wise toxicity plots.

## Stage 4: Calculate single and double mutant effects from AA PCA

This stage ('tardbpdms_aa_properties_mutant_effects') performs principal component analysis (PCA) of a curated collection of numerical indices representing various physicochemical and biochemical properties of amino acid (AA) properties. AA property feature values represent the difference between the WT and mutant PC scores.

## Stage 5: Calculate single and double mutant effects from aggregation tool predictions

This stage ('tardbpdms_agg_tools_mutant_effects') calculates aggregation / disorder algorithm feature values for single and double mutant variants (similar to stage 4).

## Stage 6: Single mutant heatmaps

This stage ('tardbpdms_single_mutant_heatmaps') produces single mutant heatmaps of toxicity effects.

## Stage 7: Human disease mutations

This stage ('tardbpdms_human_disease_mutations') tests whether human disease mutations have biased toxicity estimates.

## Stage 8: Dot plots showing explained variance of models to predict variant toxicity

This stage ('tardbpdms_toxicity_model_summary') produces plots of results from simple linear regression models to predict variant toxicity.

## Stage 9: Violin plots showing toxicity vs #introduced AAs

This stage ('tardbpdms_num_introduced_aa_violins') produces violin plots and scatterplots of toxicity (distributions) versus hydrophobicity, aromaticity, charge and aggregation propensity.

## Stage 10: Hydrophobicity of WT and toxicity hotspot line plot

This stage ('tardbpdms_wt_hydrophobicity') plots hydrophobicity score and mean toxicity effect along the length of WT TDP-43.

## Stage 11: Epistasis analysis

This stage ('tardbpdms_epistasis_analysis') performs epistasis calculations separately for each DMS library. This stage is computationally intensive (~1hour on 10 cores) and is therefore not run by default. To force re-execution of this stage set 'rerun_epistasis = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 12: Secondary structure predictions

This stage ('tardbpdms_secondary_structure_predictions') performs secondary structure predictions for each DMS library and produces combined summary plots. Secondary structure predictions are computationally intensive and are therefore not re-run by default. To force re-execution of secondary structure predictions set 'rerun_structure = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 13: Guenther structure propensities

This stage ('tardbpdms_guenther_structure_propensities') performs structure propensity calculations for PDB-structure derived
contact matrices. Structure propensity calculation are computationally intensive and are therefore not re-run by default. To force re-execution of structure propensity calculations set 'rerun_structure = T'. Additionally, the corect path to your local copy of the [DMS2structure](https://github.com/lehner-lab/DMS2structure) repository must be specified with 'DMS2structure_path = MY_LOCAL_PATH'.

## Stage 14: PWI heatmaps

This stage ('tardbpdms_PWI_heatmaps') plots pair-wise interaction (PWI) score heatmaps.


