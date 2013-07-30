MCMCDataSet<-
function(data,output.DataGeneSets,list.phenotype.ids)
		{
			cuantos.phe=length(list.phenotype.ids)	
			matriz.result=NULL
			for(i in 1:cuantos.phe)
			{
				means.dts.grps = lapply(1:length(output.DataGeneSets), function(x) GrpMean(output.DataGeneSets[x],data,list.phenotype.ids[[i]]))
				y.mu.genes = do.call(rbind, lapply(means.dts.grps, unlist))
                y.mu.genes = apply(y.mu.genes,1,function(x){mean(x,na.rm=T)})
			    matriz.result=cbind(matriz.result,y.mu.genes)
			}
			matrices=list(matriz.result)
			names(matrices)=c("y.mu")
			return(matrices)
		}