\name{GeneGrps2AffyGrpsFun}
\alias{GeneGrps2AffyGrpsFun}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that translates GO hierarchical tree of gene symbols to a groups of Affymetrix identifiers of interest.
}
\description{
This function converts GO hierarchical trees of gene symbols to groups of their associated Affymetrix identifiers.
}
\usage{
GeneGrps2AffyGrpsFun(GO.group, level, gene.u.name, affy.max.indexes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GO.group}{
The GO ontology that is of interest to consider in the analysis. ("BP","MF","CC")
}
  \item{level}{
Depth of the GO hierarchical tree of the ontology chosen. Level can only be an integer between 3 and 10.
}
  \item{gene.u.name}{
An \code{\link[GeneGroupAnalysis]{ArrayInfoFun}} item containing in a vector the unique gene symbols considered on the data set of interest.
}
  \item{affy.max.indexes}{
A \code{\link[GeneGroupAnalysis]{GeneMaxVarFun}} item consisting of the vector with the Affymetrix identifiers associated to the genes that had the maximum variability across phenotypes.
}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{index.GO.grps}{List of row number groups for their associated gene symbols.}
 \item{gene.GO.grps}{List of  gene symbol groups for their associated GO functional class.}
 \item{GO.grps}{List of GO functional classes and their associated children functional classes.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate and John Quackenbush.
}
%\note{
%
%}
\seealso{
 See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function and the help of the function \code{\link[GeneGroupAnalysis]{GibbsAllFun.ts}} for a detailed example of its use.
}
\examples{
#- For an example on the use of this function go to:
#- GibbsAllFun.ts

}
%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
