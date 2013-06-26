\name{GibbsFun.cs}
\alias{GibbsFun.cs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Gibb's sampler for the cross-sectional data design.
}
\description{
This function obtains the posterior samples for the parameters of the linear model on the mean of the gene functional classes under a cross-sectional data design.
}
\usage{
GibbsFun.cs(y.mu.a, y.mu.b, grp.sz, beta.mat, alfa.mat, sgm.y.a, sgm.y.b, sgm.alfa, rho, pi.i, mm, aa, bb, aa.pi, apriori.diff.exp, nsim, burn.in, often, prob.cut.off)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mu.a}{
A \code{\link[GeneGroupAnalysis]{MCMCData.cs}} item with the preprocessed data matrix for phenotype a.
}
  \item{y.mu.b}{
A \code{\link[GeneGroupAnalysis]{MCMCData.cs}} item with the preprocessed data matrix for phenotype b.
}
  \item{grp.sz}{
A \code{\link[GeneGroupAnalysis]{SizeGOAffyGrps}} item containing a vector with the set sizes of the gene groups with the size desired.
}
  \item{beta.mat}{
Matrix for the posterior samples for the parameter Beta.
}
  \item{alfa.mat}{
Matrix for the posterior samples for the parameter Alpha.
}
  \item{sgm.y.a}{
Matrix for the posterior samples of the variance parameter of the distribution assumed on the data for phenotype a.
}
  \item{sgm.y.b}{
Matrix for the posterior samples of the variance parameter of the distribution assumed on the data for phenotype b.
}
  \item{sgm.alfa}{
Vector for the posterior samples of the variance of the prior distribution for the parameter Alpha.
}
  \item{rho}{
Vector for the posterior samples for the parameter Rho.
}
  \item{pi.i}{
Matrix for the posterior samples for the parameter Pi.
}
  \item{mm}{
A real number for the mean of the Beta distribution for the non-zero part for the parameter Pi.
}
  \item{aa}{
A real number for the shape parameter of the inverse gamma distribution of the variance parameter of the distribution assumed on the data.
}
  \item{bb}{
A real number for the scale parameter of the inverse gamma distribution of the variance parameter of the distribution assumed on the data.
}
  \item{aa.pi}{
A real number for the precision parameter of the Beta distribution for the non-zero part for the parameter Pi.
}
  \item{apriori.diff.exp}{
An integer for the expected number of apriori differentially expressed gene groups.
}
  \item{nsim}{
An integer defining the total number of iterations for the Gibb's sampler.
}
  \item{burn.in}{
An integer defining the number of iterations that define the burn-in period for the posterior samples generated.
}
  \item{often}{
An integer defining the frequency of the iterations to print the number of differentially expressed gene functional classes up to those iterations.
}
 \item{prob.cut.off}{
A real number defining the cut off on the probability of functional gene classes differentially expressed across phenotypes. The number of functional gene classes that pass this threshold will be printed put during the iteration procedure of the Gibb's sampler.
 }
}
\details{
This function calculates the posterior samples for the parameters of interest for the cross-sectional data design, where the objective is to calculate the functional gene classes that are differentially expressed across phenotypes.
During the iteration process, this function prints out the number of differentially expressed gene functional classes with the desired probability, defined by prob.cut.off. 
}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
 \item{Beta}{Matrix with the posterior samples for the parameter beta.}
 \item{Sgma.Y.A}{Matrix with the posterior samples for the variance parameter of the distribution assumed on the data for phenotype a.}
 \item{Alfa}{Matrix with the posterior samples for the parameter alpha.} 
 \item{Sgma.Y.B}{Matrix with the posterior samples for the variance parameter of the distribution assumed on the data for phenotype b.}
 \item{Pi}{Vector with the posterior samples for the parameter pi.}
 \item{Sgm.alfa}{Vector with the posterior samples for the variance of the prior distribution for the parameter alpha.} 
 \item{Rho}{Vector with the posterior samples for the parameter rho.}
%% ...
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
A detailed example on the use of this function is provided in the \code{GeneGroupAnalysis} Vignette. 
}
\examples{
#- vignette("GeneGroupAnalysis")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
%%