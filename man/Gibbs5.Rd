\name{Gibbs5}
\alias{Gibbs5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function obtains the MCMC chains for the parameters of interest that will form their posterior distribution.
}
\description{
This function provides the MCMC chains for the parameters of interest that will form their posterior distribution. This function is to obtain the gene sets that are differentially expressed among five phenotypes of interest, taking into account one as baseline.
}
\usage{
Gibbs5(noRow,noCol,iter,GrpSzs,YMu,L0,V0,L0A,V0A,MM,AAPi,ApriDiffExp,result1,result2,result3,result4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{noRow}{
Number of row of the dataset
}
  \item{noCol}{
Total number of subjects considered.
}
  \item{iter}{
Number of iterations for the Gibbs sampler.
}
  \item{GrpSzs}{
Vector with the sizes of the gene sets considered. Output from the function \code{\link[BAGS]{DataGeneSets}}.
}
  \item{YMu}{
Output y.mu from the \code{\link[BAGS]{MCMCDataSet}}
}
  \item{L0}{
Vector with the prior parameters.
}
  \item{V0}{
Vector with the prior parameters.
}
  \item{L0A}{
Vector with the prior parameters.
}
  \item{V0A}{
Vector with the prior parameters.
}
  \item{MM}{
Parameter of the prior.
}
  \item{AAPi}{
Parameter of the prior.
}
  \item{ApriDiffExp}{
Number of differentially expressed gene sets apriori
}
  \item{result1}{
Matrix for the MCMC chains for the parameter that identifies the difference in gene set expression from phenotype 1 in comparison with the phenotype chosen as baseline. The rows are for the gene sets and the columns for the number of iterations.
}
 \item{result2}{
Matrix for the MCMC chains for the parameter that identifies the difference in gene set expression from phenotype 2 in comparison with the phenotype chosen as baseline. The rows are for the gene sets and the columns for the number of iterations.
}
 \item{result3}{
Matrix for the MCMC chains for the parameter that identifies the difference in gene set expression from phenotype 3 in comparison with the phenotype chosen as baseline. The rows are for the gene sets and the columns for the number of iterations.
}
 \item{result4}{
Matrix for the MCMC chains for the parameter that identifies the difference in gene set expression from phenotype 4 in comparison with the phenotype chosen as baseline. The rows are for the gene sets and the columns for the number of iterations.
}
}
\details{
This function provides the MCMC chains for the estimation of the posterior distribution for the parameters of interest for each gene set.
}
\value{
This function returns a list with four items
\item{alfa.1}{A list with the MCMC chains for the estimation of the posterior distribution for the parameter associated with the comparison of phenotype 1 with respect to the reference phenotype.}
\item{alfa.2}{A list with the MCMC chains for the estimation of the posterior distribution for the parameter associated with the comparison of phenotype 2 with respect to the reference phenotype.}
\item{alfa.3}{A list with the MCMC chains for the estimation of the posterior distribution for the parameter associated with the comparison of phenotype 3 with respect to the reference phenotype.}
\item{alfa.4}{A list with the MCMC chains for the estimation of the posterior distribution for the parameter associated with the comparison of phenotype 4 with respect to the reference phenotype.}
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
This function can also be used when the gene expression data has a time series experimental design. In this case, there will be five time points on the time course sampling. The assumption is that measurements between time points are independent. This assumption is reasonable when there is irregular and sparse time course sampling.
}
\examples{
# Similar to the example on Gibbs2, but in this case there are five different phenotypes of interest.  The user has to define which if the three is the reference group in order to obtain the gene groups that are differentially expressed.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
