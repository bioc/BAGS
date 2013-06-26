SgmYCFun.ts <-
function(num.grps,num.time.pnts,num.pcnts,grp.sz,y.mu.a,y.mu.b,beta.mat,alfa.mat,L0,v0,indexes,result)
            {
output=.C("conditionalSg",as.integer(num.grps),as.integer(num.time.pnts),as.integer(num.pcnts),as.integer(indexes),as.double(grp.sz),as.double(y.mu.a),as.double(y.mu.b),as.double(beta.mat),as.double(alfa.mat),as.double(L0),as.double(v0),as.double(result))
     
      Resultado=matrix(output[[12]],num.grps*num.time.pnts,num.time.pnts)
      return(Resultado)
     }

