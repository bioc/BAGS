\name{MCMCData.cs}
\alias{MCMCData.cs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that computes the required data set objects for the Gibb's sampler under a cross-sectional data design.
}
\description{
This function provides several data objects that are required in the Gibb's sampler iteration procedure. The first data objects are a transformation of the original data set into a data matrix with mean gene expression measurements where the rows of the matrix correspond to functional gene classes and the columns to the original samples provided. Finally this function provides a list of GO (Gene Ontology) processes identifiers that meet with the size requirements.
}
\usage{
MCMCData.cs(wrk.grps, data.grps, GO.grps, data, phenotype.a, phenotype.b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wrk.grps}{
A \code{\link[GeneGroupAnalysis]{SizeGOAffyGrps}} item containing the list of row numbers associated to each gene functional class.
}
  \item{data.grps}{
A \code{\link[GeneGroupAnalysis]{GeneGrps2AffyGrpsFun}} item containing the list of row numbers of the respecitve gene symbols for each gene functional class.
}
  \item{GO.grps}{
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
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{y.mu.a}{Data matrix with observations for the functional gene classes for phenotype a.}
 \item{y.mu.b}{Data matrix with observations for the functional gene classes for phenotype b.}
 \item{proc.GO}{A vector with the GO identifiers of the functional gene classes with the minimum number of elements desired.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate and John Quackenbush.
}
%%\note{
%%%%
%%}

\seealso{
A detailed example on the use of this function is provided in the \code{GeneGroupAnalysis} Vignette.  
}
\examples{
#- vignette("GeneGroupAnalysis")
}
%%% Add one or more standard keywords, see file 'KEYWORDS' in the
%%% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
