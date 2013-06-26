BetaCFun.ts <-
function(num.grps,num.time.pnts,num.pcnts,grp.sz,y.mu.a,y.mu.b,alfa.mat,sgm.y,indexes,result)
     {
output=.C("conditionalBeta",as.integer(num.grps),as.integer(num.time.pnts),as.integer(num.pcnts),as.integer(indexes),as.double(grp.sz),as.double(y.mu.a),as.double(y.mu.b),as.double(alfa.mat),as.double(sgm.y),as.double(result))
     
      Resultado=matrix(output[[10]],num.grps,num.time.pnts)
      return(Resultado)
     }

