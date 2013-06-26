AlfaFun.cs <-
function(y.mu,grp.sz,beta.mat,sgm.y,sgm.alfa,rho,mm)
            {
             results.list=lapply(1:length(grp.sz),function(x) Alfa.Aux.cs(grp.sz[x],y.mu[x,],beta.mat
[x],sgm.alfa,sgm.y[x],rho,mm))
             results=unlist(results.list)
             return(results)
            }

