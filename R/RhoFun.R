RhoFun <-
function(pi.i,ss,rr)
          {	
           I=length(pi.i)
           
           aux.1=length(which(pi.i!=0))
           shape.1=ss*rr+aux.1
           shape.2=ss*(1-rr)+(I-aux.1)
           results=rbeta(1,shape.1,shape.2)
           return(results)	
          }

