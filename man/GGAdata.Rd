\name{GGAdata}
\alias{GGAdata}
\docType{data}
\title{List containing a collection of gene symbols with their associated GO term and GO hierarchical trees for each ontology}
\description{
This dataset contains list containing a matrix with the gene symbols associated with their respective GO term, based on the  from the MSigDB at the Broad Institute. It also contains lists of GO hierarchical trees based on the ontologies and the MSigDB from the Broad Institute.}
\usage{data(GGAdata)}
\format{
  A matrix containing the following information by rows:
  \itemize{
    \item \code{GGAdata$MSigDBAnnotation[1:20,]}: The first 20 rows of the data base from MSigDD. It displays the gene symbols, GO ids, Ontology for who the GO id belongs (C-CC, F-MF, P-BP) and a description of the gene function.
    \item \code{GGAdata$CClevels[[3]]}: The CC hierarchical tree of depth 3.
    \item \code{GGAdata$MFlevels[[3]]}: The MF hierarchical tree of depth 3.
    \item \code{GGAdata$BPlevels[[3]]}: The BP hierarchical tree of depth 3.
  }
}
\details{
 This dataset contains list containing a matrix with the gene symbols associated with their respective GO term, based on the  from the MSigDB at the Broad Institute. As well it contains lists of GO hierarchical trees based on the ontologies and the MSigDB from the Broad Institute. This data set enables the construction of the gene groups of interest for the analysis of the methodology proposed.
 } 
%\itemize{
%}
\source{
\url{http://www.broadinstitute.org/gsea/msigdb/collections.jsp#C5}
}
\references{
Quiroz-Zarate A and  Quackenbush J (2012). "Manuscript in preparation".
}
\examples{
#- load the dataset
data(GGAdata,package="GeneGroupAnalysis")
#- show the first 20 rows and columns.
GGAdata$MSigDBAnnotation[1:20,]
}
\keyword{datasets}