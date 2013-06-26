SgmYFun.cs <-
function(grp.sz,y.mu,param.mat,L0,v0)
            {
             results.list=lapply(1:length(grp.sz),function(x) SgmY.Aux.cs(grp.sz[x],y.mu[x,],param.mat[x],L0,v0))	
             results=unlist(results.list)
             results
            }

