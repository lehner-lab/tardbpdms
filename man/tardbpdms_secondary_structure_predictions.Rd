% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tardbpdms_secondary_structure_predictions.R
\name{tardbpdms_secondary_structure_predictions}
\alias{tardbpdms_secondary_structure_predictions}
\title{tardbpdms_secondary_structure_predictions}
\usage{
tardbpdms_secondary_structure_predictions(toxicity_dt, outpath, miscpath,
  DMS2structure_path, colour_scheme, execute = TRUE,
  rerun_structure = F)
}
\arguments{
\item{toxicity_dt}{data.table with single and double mutant toxicity values (required)}

\item{outpath}{output path for plots and saved objects (required)}

\item{miscpath}{path to misc scripts and data directory (required)}

\item{DMS2structure_path}{Path to DMS2structure repository (required)}

\item{colour_scheme}{colour scheme file (required)}

\item{execute}{whether or not to execute the analysis (default: TRUE)}

\item{rerun_structure}{re-run structure analysis? (default:F)}
}
\value{
Nothing
}
\description{
Secondary structure predictions.
}
