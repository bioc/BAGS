Alfa.Aux.cs <-
function(grp.sz.i,y.mu,param,sgm.alfa.i,sgm.y.i,rho,mm)
                        {
                         C.SI.diff=(sgm.y.i/grp.sz.i)+sgm.alfa.i
                         C.NO.diff=sgm.y.i/grp.sz.i
                
                         Mu=y.mu-param
                         suma=dnorm(Mu,0,sqrt(C.SI.diff),log=TRUE)
                         zuma=dnorm(Mu,0,sqrt(C.NO.diff),log=TRUE)
                         rzn.log=sum(suma)-sum(zuma)+log((rho*mm)/(1-rho*mm))
                         rzn=exp(rzn.log)
                         prob=rzn/(1+rzn)
                         prob[rzn.log>600]=1
                         prob[rzn.log<(-600)]=0
                         
                         C.est=(length(y.mu)*grp.sz.i/sgm.y.i)+(1/sgm.alfa.i)
                         C.est=(1/C.est)
                         Mu2=(grp.sz.i*C.est/sgm.y.i)*sum(Mu)
                         alea=runif(1)
                         results=0
                         if(alea<prob)
                            results=rnorm(1,Mu2,sqrt(C.est))	
                         results	
                        }

