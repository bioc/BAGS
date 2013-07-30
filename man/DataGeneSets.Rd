\name{DataGeneSets}
\alias{DataGeneSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that filters the gene sets to work with the desired size.
}
\description{
This function provides the gene sets with the desired size. It provides two lists of gene sets, one with the gene identifiers of interest, and the other with the position of the gene identifiers with respect to the dataset. Finally it provides a list of the sizes of all the gene sets considered.
}
\usage{
DataGeneSets(output.ReadGMT, data.gene.symbols, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output.ReadGMT}{
Output of the function \code{\link[BAGS]{ReadGMT}}.
}
  \item{data.gene.symbols}{
Vector with the gene identifiers associated to the dataset interest. These gene identifiers have to the same as the ones on the .gmt file of interest.
}
  \item{size}{
Integer with the minimum number of genes in a gene set.
}
}
\details{
This function constructs the gene sets that are going to be considered in the analysis based on the desired size.
}
\value{
This function returns a list with three items
\item{DataGeneSetsIds}{A list of gene sets with the positions of the gene identifiers with respect to the dataset of interest.}
\item{DataGeneSetsNms}{A list of gene sets of the gene identifiers.}
\item{Size}{A vector with the size of the gene sets}
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

data(vdx,package="breastCancerVDX")
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


#===== In the following R dataset it is stored the .gmt file associated to the MF from GO.
#===== So "reading the GMT" is the only step that we skip. But an example is provided on the
#===== help file associated to the function "ReadGMT".
data(AnnotationMFGO,package="BAGS")

data.gene.grps=DataGeneSets(AnnotationMFGO,gene.nams,10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
