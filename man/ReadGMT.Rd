\name{ReadGMT}
\alias{ReadGMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that extracts the gene set definition from a .gmt file.
}
\description{
This function reads the gene set definitions provided by a .gmt file (MSigDB annotation files) an stores the information into a list.
}
\usage{
ReadGMT(path)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
The path where the .gmt file is stored..
}
}
\details{
This function reads the gene set definitions provided by a .gmt file (MSigDB annotation files) an stores the information into a list. The .gmt file needs to have the name of the gene set in column 1. A description/properties of the gene set in column 2. And finally the gene identifiers associated to the gene set from column 3 and on. These gene identifiers can be Ensembl, Affymetrix, Gene Symbols, etc...
}
\value{
This function returns a list of the same length as gene sets in the .gmt file provided. The entries of the lists are the associated gene identifier associated to their respective gene set
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
\author{
A. Quiroz-Zarate
\email{aquiroz@jimmy.harvard.edu}
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See the \code{BAGS} Vignette for examples on how to use function \code{\link[BAGS]{Gibbs2}}.
}
\examples{
# An example on the use of this function:
# gene.sets.DB=ReadGMT("/Users/Bioinformatics/Projects/c5.mf.v3.0.symbols.gmt")
#
# #Where the path is where the gmt file from the MSigDB is stored.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gmt}
\keyword{MSigDB}% __ONLY ONE__ keyword per line
