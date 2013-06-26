\name{MCMCData.ts}
\alias{MCMCData.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that computes the required data set objects for the Gibb's sampler under a time series data design.
}
\description{
This function provides several data objects that are required in the Gibb's sampler iteration procedure. The first data objects are a transformation of the original data set into a data matrix with mean gene expression measurements where the rows of the matrix correspond to functional gene classes and the columns to the original samples provided. It also provides the covariance matrix required for the posterior sample calculation of the covariance of the data. Finally this function provides a list of GO (Gene Ontology) processes identifiers that meet with the size requirements.
}
\usage{
MCMCData.ts(wrk.grps, data.grps, GO.proc, data, phenotype.a, phenotype.b, indexes, num.time.pnts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wrk.grps}{
A \code{\link[GeneGroupAnalysis]{SizeGOAffyGrps}} item containing the list of row numbers associated to each gene functional class.
}
  \item{data.grps}{
A \code{\link[GeneGroupAnalysis]{GeneGrps2AffyGrpsFun}} item containing the list of row numbers of the respective gene symbols for each gene functional class.
}
  \item{GO.proc}{
A \code{\link[GeneGroupAnalysis]{GeneGrps2AffyGrpsFun}} item containing the list of GO gene functional classes.
}
  \item{data}{
Data set of interest in the form of a matrix: columns for patients, rows for Affymetrix identifiers.
}
  \item{phenotype.a}{
Vector of indices with the column identifiers of phenotype a.
}
  \item{phenotype.b}{
Vector of indices with the column identifiers of phenotype b.
}
  \item{indexes}{
Vector with the order in which the subjects are recorded across the time course in the data set.
}
  \item{num.time.pnts}{
Number of time points considered.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{y.mu.a}{Data matrix with observations for the functional gene classes for phenotype a.}
 \item{y.mu.b}{Data matrix with observations for the functional gene classes for phenotype b.}
 \item{lambda}{A priori matrix estimate for the covariance matrix to be used for the posterior sample estimates for the variance of the distribution assumed on the data.}
 \item{proc.GO}{A vector with the GO identifiers of the functional gene classes with the minimum number of elements desired.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate and John Quackenbush.
}
\seealso{
 See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function and the help of the function \code{\link[GeneGroupAnalysis]{GibbsAllFun.ts}} for a detailed example of its use.
}
\examples{
#- For an example on the use of this function go to:
#- GibbsAllFun.ts

}