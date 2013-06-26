SgmY.Aux.cs <-
function(grp.sz.i,y.mu,param,L0,v0)
              {
               M=sum((y.mu-param)^2)          
               scala=(grp.sz.i*M/2)+L0
               forma=(length(y.mu)/2)+v0
               y=rgamma(1,shape=forma,rate=scala)                 
               z=1/y
               return(z)	
              }

