BetaFun.cs <-
function(y.mu.a,y.mu.b,grp.sz,alfa.vec,sgm.y.a,sgm.y.b)
                {
                 results.list=lapply(1:length(grp.sz),function(x) Beta.Aux.cs(grp.sz[x],y.mu.a[x,],y.mu.b[x,],alfa.vec[x],sgm.y.a[x],sgm.y.b[x]))
                 results=unlist(results.list)
                 return(results)
                }

