MCMCData.cs <-
function(wrk.grps,data.grps,GO.grps,data,phenotype.a,phenotype.b)
                {	
                 means.dts.grps.1=lapply(1:length(wrk.grps),function(x) GrpMeanDataFun(data.grps[wrk.grps[x]],data,phenotype.a))
                 means.dts.grps.2=lapply(1:length(wrk.grps),function(x) GrpMeanDataFun(data.grps[wrk.grps[x]],data,phenotype.b))
                 y.mu.genes.1=do.call(rbind,lapply(means.dts.grps.1,unlist))
                 y.mu.genes.2=do.call(rbind,lapply(means.dts.grps.2,unlist))
                 p.M.t=GO.grps[wrk.grps]
                 resultds=list(y.mu.genes.1,y.mu.genes.2,p.M.t)
                 names(resultds)=c("y.mu.a","y.mu.b","proc.GO")
                 return(resultds) 
                }

