SgmAlfaFun.cs <-
function(alfa.vec,aa,bb)
                    {
                     I=length(alfa.vec)
                     forma=(I/2)+aa
                     scala=(sum(alfa.vec^2)/2)+bb	
                     y=rgamma(1,shape=forma,rate=scala)
                     resultado=1/y
                     return(resultado)
                    }

