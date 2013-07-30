\name{MCMCDataSet}
\alias{MCMCDataSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that transform the dataset into the required format for the Gibbs sampler.
}
\description{
This function makes the necessary transformation of the dataset in order for the Gibbs sampler to perform the iterations. This transformation is based on the number of phenotypes of interest considered.
}
\usage{
MCMCDataSet(data,output.DataGeneSets,list.phenotype.ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The dataset of interest. The data has as rows the gene identifiers of interest (same as the .gmt file) and as columns the samples considered.
}
  \item{output.DataGeneSets}{
List of gene sets with the positions of the gene identifiers with respect to the dataset of interest. This list is part of the output from \code{\link[BAGS]{DataGeneSets}}.
}
  \item{list.phenotype.ids}{
A list that has as elements the vectors with the column positions of the phenotypes considered in the analysis.
}
}
\details{
This function constructs the gene sets that are going to be considered in the analysis based on the desired size.
}
\value{
This function returns a list:
\item{y.mu}{A matrix with the means across samples and genes for each gene set. The rows are gene sets and columns are the different phenotypes considered.}
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
See the \code{BAGS} Vignette for examples on how to use this function and the help of the function \code{\link[BAGS]{Gibbs5}} for a detailed example of its use.
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

data.gene.grps=DataGeneSets(AnnotationMFGO,gene.nams,10)

phntp.list=list(er.pos,er.neg)
data.mcmc=MCMCDataSet(gene.expr,data.gene.grps$DataGeneSetsIds,phntp.list)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
