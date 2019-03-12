# Overview

Welcome to the GitHub repository for the following publication: "The mutational landscape of a prion-like domain" (Bolognesi B & Faure AJ et al., 2019)

Here you'll find an R package with all scripts to reproduce the figures and results from the computational analyses described in the paper.

# Required Software

To run the tardbpdms pipeline you will need the following software and associated packages:

* **[R](https://www.r-project.org/) >=v3.5.2** (data.table, ggplot2, parallel, plyr, reshape2, lemon, Biostrings, GGally, ppcor)

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

All pre-processed data and required miscellaneous files first need to be downloaded from [here](https://www.dropbox.com/s/nbuvtfler395rcd/tardbpdms_misc.zip?dl=0).


