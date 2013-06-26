\name{SizeGOAffyGrps}
\alias{SizeGOAffyGrps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute gene functional classes with a desired size.
}
\description{
This obtains the gene functional classes with the minimum desired size as well as providing the size of such functional classes. 
}
\usage{
SizeGOAffyGrps(affy.GO.grps, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{affy.GO.grps}{
A \code{\link[GeneGroupAnalysis]{GeneGrps2AffyGrpsFun}} item with a list of groups with the row number of their respective Affymetrix identifiers of the gene functional classes.
}
  \item{size}{
An integer defining the minimum gene functional size desired.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{groups}{A list of  groups with the row numbers of the respective Affymetrix identifiers of the the functional gene classes, with the group size desired.}
 \item{group.size}{Vector with the set sizes of the groups with the size desired.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate and John Quackenbush.
}
%%\note{
%%  ~~further notes~~
%%}
\seealso{
 See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function and the help of the function \code{\link[GeneGroupAnalysis]{GibbsAllFun.ts}} for a detailed example of its use.
}
\examples{
#- For an example on the use of this function go to:
#- GibbsAllFun.ts

}
%%% Add one or more standard keywords, see file 'KEYWORDS' in the
%%% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line