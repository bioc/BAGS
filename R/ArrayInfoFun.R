ArrayInfoFun <-
function(data,array.symbols)
                   {
					array.list=AnnotationDbi::as.list(array.symbols)
                   	data.order=match(rownames(data),names(array.list),nomatch=0)
                   	which.no.na=which(unlist(array.list)[data.order]!="NA")
                   	genes.name=unlist(array.list)[data.order][which.no.na]
                    probes.name=names(genes.name)
                   	genes.name.unique=unique(genes.name)
                   	probes.ids=data.order[which.no.na]
                    results=list(genes.name,probes.name,probes.ids,genes.name.unique)
                    names(results)=c("genes.name","probes.name","probes.ids","genes.name.unique")
                    return(results)
                   }

