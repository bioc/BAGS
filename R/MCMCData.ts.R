MCMCData.ts <-
function(wrk.grps,data.grps,GO.proc,data,phenotype.a,phenotype.b,indexes,num.time.pnts)
                {
                 num.pcnts=length(indexes)		
                 mean.dts.grps.1=lapply(1:length(wrk.grps),function(x) GrpMeanDataFun(data.grps[wrk.grps[x]],data,phenotype.a))
                 mean.dts.grps.2=lapply(1:length(wrk.grps),function(x) GrpMeanDataFun(data.grps[wrk.grps[x]],data,phenotype.b))
                 y.mu.genes.1=do.call(rbind,lapply(mean.dts.grps.1,unlist))
                 y.mu.genes.2=do.call(rbind,lapply(mean.dts.grps.2,unlist))
                 
                 SgmAprGrps.ts=lapply(1:length(wrk.grps),function(x) SgmAprGrps.Aux.ts(y.mu.genes.1[x,],y.mu.genes.2[x,],indexes,num.time.pnts,num.pcnts))
                 SgmAprGrps.ts=MeanMatGrps.Aux.ts(SgmAprGrps.ts,num.time.pnts)
                 p.M.t=GO.proc[wrk.grps]
                 
                 resultds=list(y.mu.genes.1,y.mu.genes.2,SgmAprGrps.ts,p.M.t)
                 names(resultds)=c("y.mu.a","y.mu.b","lambda","proc.GO")
                 return(resultds) 
                }

