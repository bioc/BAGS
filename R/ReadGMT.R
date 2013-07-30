ReadGMT<-
function(path)
        {
        	#=="input": path of the .gmt file of interest
        	#=="output": list of the gene groups of interest
        	gene.sets.lines=readLines(path)
        	gene.sets.list=lapply(gene.sets.lines,function(x) {unlist(strsplit(x, "\t"))})
        	gene.ids.list.names=unlist(lapply(gene.sets.list,function(x){aux=x[1]}))
	        gene.ids.list.grps=lapply(gene.sets.list,function(x){aux=x[-c(1:2)]})              
	        names(gene.ids.list.grps)=gene.ids.list.names
	        return(gene.ids.list.grps)
        }