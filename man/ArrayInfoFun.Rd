\name{ArrayInfoFun}
\alias{ArrayInfoFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that extracts the required information from the data set.
}
\description{
This function captures the Affymetrix identifiers and the gene symbols associated to the data set of interest in order to prepare the data set for the Gibb's sampler execution.
}
\usage{
ArrayInfoFun(data, array.symbols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Data set of interest in the form of a matrix: columns for patients, rows for Affymetrix identifiers.
}
  \item{array.symbols}{
Mapping between Affymetrix identifiers and gene symbols for the associated Affymetrix platform of the data set of interest.
}
}
\details{
This function captures the necessary information from the data set of interest.
}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
\item{genes.name}{A vector of gene symbols off all the Affymetrix identifiers (rows of the data set) that have a gene symbol.}
\item{probes.name}{A vector of Affymetrix identifiers (rows of the data set).}
\item{probes.ids}{A vector containing the row number of their respective Affymetrix identifier (rows of the data set).}
\item{genes.name.unique}{A vector of unique gene symbols off all the Affymetrix identifiers (rows of the data set) in the data set.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate and John Quackenbush.
}
%%\note{
%% See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function. 
%%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function and the help of the function \code{\link[GeneGroupAnalysis]{GibbsAllFun.ts}} for a detailed example of its use.
}
\examples{
#- For an example on the use of this function go to:
#- GibbsAllFun.ts
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
