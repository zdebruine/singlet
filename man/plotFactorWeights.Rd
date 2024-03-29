% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotFactorWeights.R
\name{plotFactorWeights}
\alias{plotFactorWeights}
\title{convenience function to map one or more factors along a genome using igvR}
\usage{
plotFactorWeights(object, gr, factors = 1:3, plot = FALSE)
}
\arguments{
\item{object}{an nmf object or something with a @w weights matrix}

\item{gr}{a GRanges object with coordinates for the features}

\item{factors}{which factors to plot weights for (default: 1, 2, 3)}

\item{plot}{use igvR to plot the factors? (TRUE, if igvR detected)}
}
\value{
the GRanges gr, but with factor weights added as mcols
}
\description{
convenience function to map one or more factors along a genome using igvR
}
\details{
This function presumes a GRanges object will be supplied, which in turn
 presumes that the GenomicRanges package is installed from Bioconductor. 
 Further, if plot == TRUE, the igvR package is presumed to be installed. 
 If either of these presumptions are false, or if factor weights cannot
 be mapped to identifiers in the GRanges, this function will fail.
}
