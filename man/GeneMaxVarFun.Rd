\name{GeneMaxVarFun}
\alias{GeneMaxVarFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Affymetrix identifiers with the maximum variability between phenotypes.
}
\description{
Some Affymetrix identifiers map to the same gene symbol. This implies that there exist more than one row number from the data set associated to a gene symbol.
This function computes the respective row number associated to the Affymetrix identifier the presents the maximum variability across phenotypes for every group of Affymetrix identifiers mapping to gene symbols.
}
\usage{
GeneMaxVarFun(ArrayInfoFun.ans, data, phenotype.a, phenotype.b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ArrayInfoFun.ans}{
A \code{\link[GeneGroupAnalysis]{ArrayInfoFun}} item of the data set of interest.
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
\details{
This function captures the row numbers associated to the Affymetrix identifiers with the maximum variability between phenotypes for every group of Affymetrix identifiers mapping to gene symbols.
}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{IndexMaxVar}{A vector with the row number ids of the respective Affymetrix identifiers with the maximum variability across phenotypes.}
 \item{AffyIdsMaxVar}{A vector with the Affymetrix identifiers associated to the "indexes.good"}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
