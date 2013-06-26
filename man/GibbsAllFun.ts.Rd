\name{GibbsAllFun.ts}
\alias{GibbsAllFun.ts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Gibb's sampler for the series data design.
}
\description{
This function obtains the posterior samples for the parameters of the linear model on the mean of the gene functional classes under a time series data design. These posterior samples are for the parameters for the model that identifies gene functional classes that present a parallel and non-parallel time course trajectories between phenotypes.
}
\usage{
GibbsAllFun.ts(y.mu.a, y.mu.b, grp.sz, beta.mat, alfa.mat, sgm.alfa, rho, pi.i, mm, aa.pi, lmbd, df.lmbd, lmbd.alf, indexes, num.time.pnts, num.pcnts, apriori.diff.exp, nsim, burn.in, often, prob.cut.off)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y.mu.a}{
A \code{\link[GeneGroupAnalysis]{MCMCData.ts}} item with the preprocessed data matrix for phenotype a.
}
  \item{y.mu.b}{
A \code{\link[GeneGroupAnalysis]{MCMCData.ts}} item with the preprocessed data matrix for phenotype b.
}
  \item{grp.sz}{
A \code{\link[GeneGroupAnalysis]{MCMCData.ts}} item with the preprocessed data matrix for phenotype b.
}
  \item{beta.mat}{
Matrix for the posterior samples for the parameter Beta.
}
  \item{alfa.mat}{
Matrix for the posterior samples for the parameter Alpha.
}
  \item{sgm.alfa}{
Matrix for the posterior samples of the variance of the prior distribution for the parameter Alpha.
}
  \item{rho}{
Vector for the posterior samples for the parameter Rho.
}
  \item{pi.i}{
Vector for the posteriors samples for the parameter Pi.
}
  \item{mm}{
A real number for the mean of the Beta distribution for the non-zero part for the parameter Pi.
}
  \item{aa.pi}{
A real number for the precision parameter of the Beta distribution for the non-zero part for the parameter Pi.
}
  \item{lmbd}{
A prior estimate for the distribution for the variance assumed on the data.
}
  \item{df.lmbd}{
An integer defining the degrees of freedom for the matrix apriori estimate for the distribution for the variance assumed on the data.
}
  \item{lmbd.alf}{
A prior estimate for the prior distribution for the variance for the parameter Alpha.
}
  \item{indexes}{
Vector of integers defining the indexes in the data that delineate the different time points for the sample identifiers (columns of the data set).
}
  \item{num.time.pnts}{
An integer defining the number of time points considered.
}
  \item{num.pcnts}{
An integer defining the number of subjects considered.
}
  \item{apriori.diff.exp}{
An integer for the expected number of apriori differentially expressed gene groups.
}
  \item{nsim}{
An integer defining the total number of iterations for the Gibb's sampler.
}
  \item{burn.in}{
An integer defining the total number of iterations for the Gibb's sampler.
}
  \item{often}{
An integer defining the frequency of the iterations to print the number of differentially expressed gene functional classes up to those iterations.
}
\item{prob.cut.off}{
A real number defining the cut off on the probability of functional gene classes differentially expressed between phenotypes across the time course. The number of functional gene classes that pass this threshold will be printed put during the iteration procedure of the Gibb's sampler.
}}
\details{
This function calculates the posterior samples for the parameters of interest for the time series data design, where the objective is to calculate the functional gene classes that are differentially expressed between phenotypes across the time course. 
During the iteration process, this function prints out the number of differentially expressed gene functional classes with the desired probability, defined by prob.cut.off.  
}
\value{
This function returns a list containing the following results:
%%  If it is a LIST, use
  \item{Beta}{Matrix with the posterior samples for the parameter beta.}
  \item{Alfa}{Matrix with the posterior samples for the parameter alpha.}
  \item{Pi}{Vector with the posterior samples for the parameter pi.}
  \item{Rho}{Vector with the posterior samples for the parameter rho.}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 See the \code{GeneGroupAnalysis} Vignette for examples on how to use this function. 
}
\examples{
library(annotate)
library(rheumaticConditionWOLLBOLD)
library(hgu133plus2.db)

data(wollbold,package="rheumaticConditionWOLLBOLD")
#----Normalized expression data set
woll.data.exp=exprs(wollbold)
#---- Phenotypic information from the data set
#pData(wollbold)
#---ER status of patients
TGF.exp=1:30
TNF.exp=31:60

#--- Checking that the columns correspond to their respective phenotype data id
all(colnames(woll.data.exp)==rownames(pData(wollbold)))
#---Checking the annotation of the data
annotation(wollbold)
array.info=ArrayInfoFun(woll.data.exp,hgu133plus2SYMBOL)
genes.max.var=GeneMaxVarFun(array.info,woll.data.exp,TGF.exp,TNF.exp)
GO2Gene.grps=GeneGrps2AffyGrpsFun("CC",3,array.info$genes.name.unique,genes.max.var$IndexMaxVar)
GO.grps=GO2Gene.grps$GO.grps

Wollbold09WrkngGrps=SizeGOAffyGrps(GO2Gene.grps$index.GO.grps,100)
indexes.1=c(1,2,3,4,5,6)
Wollbold09MCMCData=MCMCData.ts(Wollbold09WrkngGrps$groups,GO2Gene.grps$index.GO.grps,
GO2Gene.grps$GO.grps,woll.data.exp,TGF.exp,TNF.exp,indexes.1,5)

nsim = 40
burn.in = 10
Grps.apriori.diff.exp = 23
shape = 3
scale = 0.1
mm.pi = 0.75
aa.pi = 10
often = 20
cut.off = 0.7
df.lambda = 10
no.time.pnts = 5
no.pcnts = 6
v.SS.i = array(0, c(length(Wollbold09MCMCData$proc.GO), no.time.pnts, no.time.pnts,8,2))
v.SS.Al = array(0, c(no.time.pnts , no.time.pnts , 2))
v.beta.i = array(0, c(length(Wollbold09MCMCData$proc.GO), no.time.pnts, nsim))
v.alfa.i = array(0, c(length(Wollbold09MCMCData$proc.GO), no.time.pnts, nsim))
v.pi.a.i = matrix(0, length(Wollbold09MCMCData$proc.GO), nsim)
v.rho.a = rep(0, nsim)
v.beta.i[, , 1] = rnorm(length(Wollbold09MCMCData$proc.GO), 0, 1)
v.alfa.i[, , 1] = rnorm(length(Wollbold09MCMCData$proc.GO), 0, 1)
v.SS.Al[, , 1] = diag(0.1, no.time.pnts)
v.pi.a.i[, 1] = runif(length(Wollbold09MCMCData$proc.GO))
v.rho.a[1] = 0.1
L0.alfa = diag(0.005,no.time.pnts)
indexes.1 = c(0, 6, 12, 18, 24)

results.MCMC = GibbsAllFun.ts(Wollbold09MCMCData$y.mu.a, Wollbold09MCMCData$y.mu.b, Wollbold09WrkngGrps$group.size, v.beta.i, v.alfa.i, v.SS.Al, v.rho.a, v.pi.a.i, mm.pi, aa.pi, Wollbold09MCMCData$lambda, df.lambda, L0.alfa, indexes.1, no.time.pnts, no.pcnts, Grps.apriori.diff.exp, nsim, burn.in, often, cut.off)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
