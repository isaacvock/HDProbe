% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HDProbe.R
\name{VarianceTrend}
\alias{VarianceTrend}
\title{Estimate replicate variability trend}
\usage{
VarianceTrend(Filter_df, Homosked = FALSE)
}
\arguments{
\item{Filter_df}{a dataframe}

\item{Homosked}{Logical; if TRUE, homoskedasticity is assumed}
}
\value{
a list of linear model fits
}
\description{
Estimate replicate variability trend
}
