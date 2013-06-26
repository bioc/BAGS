SgmAlfaFun.ts <-
function(alfa.mat,L0,v0)
                 {
                  require(MCMCpack)	
                  
                  I=dim(alfa.mat)[1]
                  J=dim(alfa.mat)[2]	
                  unos=which(alfa.mat[,1]!=0)
                  vn=v0
                  L=L0
                  if(length(unos)>0)
                    {
                     vn=(v0+length(unos))
                     SS.Al=matrix(0,J,J)
                        for(i in 1:length(unos))
                            SS.Al=SS.Al+alfa.mat[unos[i],]%*%t(alfa.mat[unos[i],])
                     L=L0+SS.Al
                    } 
                  Ln=solve(L)
                  results=riwish(vn,Ln) 
                  results
                 }

