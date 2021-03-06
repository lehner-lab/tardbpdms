% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tardbpdms_human_disease_mutations.R
\name{tardbpdms_human_disease_mutations}
\alias{tardbpdms_human_disease_mutations}
\title{tardbpdms_human_disease_mutations}
\usage{
tardbpdms_human_disease_mutations(toxicity_dt, missense_AF_file,
  disease_mut_file, outpath, colour_scheme, execute = TRUE)
}
\arguments{
\item{toxicity_dt}{data.table with single mutant toxicity values (required)}

\item{missense_AF_file}{table of missense mutation allele frequencies (required)}

\item{disease_mut_file}{table of human disease mutations and classifications (required)}

\item{outpath}{output path for plots and saved objects (required)}

\item{colour_scheme}{colour scheme file (required)}

\item{execute}{whether or not to execute the analysis (default: TRUE)}
}
\value{
Nothing
}
\description{
Check fitness bias of human disease mutations.
}
