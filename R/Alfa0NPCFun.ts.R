Alfa0NPCFun.ts <-
function(num.grps,num.time.pnts,num.pcnts,grp.sz,y.mu,beta.mat,alfa.mat,sgm.y,indexes,result)
     {
output=.C("conditionalNPAlfa0",as.integer(num.grps),as.integer(num.time.pnts),as.integer(num.pcnts),as.integer(indexes),as.double(grp.sz),as.double(y.mu),as.double(beta.mat),as.double(alfa.mat),as.double(sgm.y),as.double(result))
     
      resultado=matrix(output[[10]],num.grps,1)
      return(resultado)
     }

