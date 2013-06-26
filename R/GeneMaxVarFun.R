GeneMaxVarFun<-
function(ArrayInfoFun.ans,data,phenotype.a,phenotype.b)
                     {  
                      genes.name.u=ArrayInfoFun.ans[[4]]
                      genes.name.a=ArrayInfoFun.ans[[1]]
                      probes.ids=ArrayInfoFun.ans[[3]]
                      
                      aux=match(genes.name.a,genes.name.u)
                      indices=0
                      for(i in 1:length(genes.name.u))
                         {
                          cuales=probes.ids[which(aux==i)]
                          if(length(cuales)>1)
                            {
                             var.r=apply(cbind(apply(data[cuales,phenotype.a],1,mean),apply(data[cuales,phenotype.b],1,mean)),1,var)
                             cuales=cuales[which.max(var.r)]
                            }
                          indices=c(indices,cuales)  	
                         }
                      indices=indices[-1]
                      affy.ids=rownames(data)[indices]
	                  results=list(indices,affy.ids)
	                  names(results)=c("IndexMaxVar","AffyIdsMaxVar")
	                  return(results)
                     }
