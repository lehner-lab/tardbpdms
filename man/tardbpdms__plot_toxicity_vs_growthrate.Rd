% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tardbpdms__plot_toxicity_vs_growthrate.R
\name{tardbpdms__plot_toxicity_vs_growthrate}
\alias{tardbpdms__plot_toxicity_vs_growthrate}
\title{tardbpdms__plot_toxicity_vs_growthrate}
\usage{
tardbpdms__plot_toxicity_vs_growthrate(input_dt, output_file,
  colour_scheme = NULL, width = 4, height = 3)
}
\arguments{
\item{input_dt}{a data.table with columns: growth_rate, growth_rate_low, growth_rate_high, toxicity, sigma, region (required)}

\item{output_file}{path to output plot (required)}

\item{colour_scheme}{colour scheme list}

\item{width}{plot width}

\item{height}{plot height}
}
\value{
Nothing
}
\description{
Plot toxicity versus growth rate.
}
