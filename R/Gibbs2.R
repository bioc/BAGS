Gibbs2<-
function(noRow,noCol,iter,GrpSzs,YMu,L0,V0,L0A,V0A,MM,AAPi,ApriDiffExp,result1)
     {
output=.C("Gibbs2",as.integer(noRow),as.double(noCol),as.integer(iter),as.double(GrpSzs),as.double(YMu),as.double(L0),as.double(V0),as.double(L0A),as.double(V0A),as.double(MM),as.double(AAPi),as.double(ApriDiffExp),as.double(result1))
      Resultado1=matrix(output[[13]],noRow,iter)
      
      Resultado=list(Resultado1)
      names(Resultado)=c("alfa.1")
      return(Resultado)
     }