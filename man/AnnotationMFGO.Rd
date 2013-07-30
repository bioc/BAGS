\name{AnnotationMFGO}
\alias{AnnotationMFGO}
\docType{data}
\title{List containing a collection of gene symbols with their associated GO term for the Molecular function ontology}
\description{
This dataset contains list containing all the GO terms for the Molecular function ontology with their associated gene symbols, base on the collections provided by the MSigDB from the Broad Institute}
\usage{data(AnnotationMFGO)}
\format{
  A list containing the information in the following way:
  \itemize{
    \item \code{AnnotationMFGO[1:3]}: The first 3 GO terms  of the Molecular Function ontology with their associated gene symbols.
  }
}
\details{
 This dataset contains list containing gene symbols associated with their respective GO term, based on the  from the MSigDB at the Broad Institute. This dataset enables the construction of the gene groups of interest for the analysis of the methodology proposed.
 } 
%\itemize{
%}
\source{
\url{http://www.broadinstitute.org/gsea/msigdb/collections.jsp#C5}
}
\references{
Quiroz-Zarate A and  Quackenbush J (2012). "Manuscript in preparation".
}
\author{
A. Quiroz-Zarate.
\email{aquiroz@jimmy.harvard.edu}
}
\examples{
#- load the dataset
data(AnnotationMFGO,package="BAGS")
#- show the first 20 rows and columns.
AnnotationMFGO[1:3]
}
\keyword{Molecular function}
\keyword{Gene Ontology}
\keyword{MSigDB}