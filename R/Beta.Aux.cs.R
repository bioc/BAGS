Beta.Aux.cs <-
function(grp.sz.i,y.mu.a,y.mu.b,param,sgm.y.a.i,sgm.y.b.i)
                        {
                         M.a=sum(y.mu.a)/sgm.y.a.i
                         M.b=sum((y.mu.b-param))/sgm.y.b.i
                         Mu2=M.a+M.b
                         C.A.est=(length(y.mu.a)*grp.sz.i/sgm.y.a.i)
                         C.B.est=(length(y.mu.b)*grp.sz.i/sgm.y.b.i)
                         C.est=C.A.est+C.B.est
                         C.est=(1/C.est)
                         Mu=grp.sz.i*C.est*Mu2
                         res=rnorm(1,Mu,sqrt(C.est))
                         return(res)
                        }

