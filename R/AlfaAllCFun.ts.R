AlfaAllCFun.ts <-
function(num.grps,num.time.pnts,num.pcnts,grp.sz,y.mu,beta.mat,sgm.y,sgm.alfa,rho,mm,indexes,result)
     {
output=.C("conditionalALLAlfa",as.integer(num.grps),as.integer(num.time.pnts),as.integer(num.pcnts),as.integer(indexes),as.double(grp.sz),as.double(y.mu),as.double(beta.mat),as.double(sgm.y),as.double(sgm.alfa),as.double(rho),as.double(mm),as.double(result))
     
      Resultado=matrix(output[[12]],num.grps,num.time.pnts)
      return(Resultado)
     }

