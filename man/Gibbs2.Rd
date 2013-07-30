\name{Gibbs2}
\alias{Gibbs2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function obtains the MCMC chains for the parameters of interest that will form their posterior distribution.
}
\description{
This function provides the MCMC chains for the parameters of interest that will form their posterior distribution. This function is to obtain the gene sets that are differentially expressed among five phenotypes of interest, taking into account one as baseline.
}
\usage{
Gibbs2(noRow,noCol,iter,GrpSzs,YMu,L0,V0,L0A,V0A,MM,AAPi,ApriDiffExp,result1)
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
Matrix for the MCMC chains for the parameter that identifies the difference in geneset expression from phenotype 1 in comparison with the phenotype chosen as baseline. The rows are for the gene sets and the columns for the number of iterations.
}
}
\details{
This function provides the MCMC chains for the estimation of the posterior distribution for the parameters of interest for each gene set.
}
\value{
This function returns a list with four items
\item{alfa.1}{A list with the MCMC chains for the estimation of the posterior distribution for the parameter associated with the comparison of phenotype 1 with respect to the phenotype chosen as baseline.}
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
See the \code{BAGS} Vignette for examples on how to use this function and the help of the function \code{\link[BAGS]{Gibbs2}} for a detailed example of its use.
This function can also be used when the gene expression data has a time series experimental design. In this case, there will be two time points on the time course sampling. The assumption is that measurements between time points are independent. This assumption is reasonable when there is irregular and sparse time course sampling.
}
\examples{
library(breastCancerVDX)
library(Biobase)

data(vdx)
gene.expr=exprs(vdx)   # Gene expression of the package
vdx.annot=fData(vdx)   # Annotation associated to the dataset
vdx.clinc=pData(vdx)   # Clinical information associated to the dataset 

# Identifying the sample identifiers associated to ER+ and ER- breast cancer
er.pos=which(vdx.clinc$er==1)
er.neg=which(vdx.clinc$er==0)

# Only keep columns 1 and 3, probeset identifiers and Gene symbols respectively
vdx.annot=vdx.annot[,c(1,3)]

all(rownames(gene.expr)==as.character(vdx.annot[,1]))  # Checking if the probeset are ordered with respect to the dataset
all(colnames(gene.expr)==as.character(vdx.clinc[,1]))  # Checking if the sample identifiers are order with respect to the dataset
rownames(gene.expr)=as.character(vdx.annot[,2])        # Changing the row identifiers to the gene identifiers of interest

#===== Because we have several measurements for a gene (multiple rows for a gene), we filter the genes
#===== Function to obtain the genes with highest variabilty among phenotypes
gene.nms.u=unique(rownames(gene.expr))
gene.nms=rownames(gene.expr)
indices=NULL
for(i in 1:length(gene.nms.u))
{
	aux=which(gene.nms==gene.nms.u[i])
	if(length(aux)>1){
		var.r = apply(cbind(apply(gene.expr[aux,er.pos],1,mean),apply(gene.expr[aux,er.neg],1,mean)),1,var)
		aux=aux[which.max(var.r)]
	}
	indices=c(indices,aux)
}
#===== Only keep the genes with most variability among the phenotypes of interest
gene.expr=gene.expr[indices,]
gene.nams=rownames(gene.expr)     # The gene symbols of interest are stored here

dim(gene.expr)

#===== In the following R dataset it is stored the .gmt file associated to the MF from GO.
#===== So "reading the GMT" is the only step that we skip. But an example is provided on the
#===== help file associated to the function "ReadGMT".
data(AnnotationMFGO,package="BAGS")

data.gene.grps=DataGeneSets(AnnotationMFGO,gene.nams,5)

phntp.list=list(er.pos,er.neg)
data.mcmc=MCMCDataSet(gene.expr,data.gene.grps$DataGeneSetsIds,phntp.list)

noRow=dim(data.mcmc$y.mu)[1]
noCol=unlist(lapply(phntp.list,length))
iter=10000
GrpSzs=data.gene.grps$Size
YMu=data.mcmc$y.mu
L0=rep(2,2)
V0=rep(3,2)
L0A=rep(3,2)
V0A=rep(3,2)
MM=0.55
AAPi=10
ApriDiffExp=floor(dim(data.mcmc$y.mu)[1]*0.03)
results=matrix(0,noRow,iter)
		
mcmc.chains=Gibbs2(noRow,noCol,iter,GrpSzs,YMu,L0,V0,L0A,V0A,MM,AAPi,ApriDiffExp,results)
		
burn.in=2000
alfa.pi=apply(mcmc.chains[[1]][,burn.in:iter],1,function(x){
	                                y=length(which(x!=0))/length(burn.in:iter);return(y)})
plot(alfa.pi,type="h",main="Probabilities of MF differentially expressed between ER status")                                   
cut.off=0.9
abline(h=cut.off,col="red")
differential.processes=names(data.gene.grps$Size)[which(alfa.pi>cut.off)]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
