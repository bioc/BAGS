PiFun <-
function(alfa.vec,rho,aa,mm)
            {
             results=rep(0,length(alfa.vec))	              
             shape.1=aa*mm
             shape.2=aa*(1-mm)	
             
             zeros=which(alfa.vec==0)
             no.zeros=which(alfa.vec!=0)
             c.0=(1-mm)*rho/(1-rho)
             prob=c.0/(1+c.0)
             alea=runif(length(zeros))
          
             results[zeros[alea<prob]]=rbeta(length(zeros[alea<prob]),shape.1,(shape.2+1))
             results[no.zeros]=rbeta(length(no.zeros),(shape.1+1),shape.2) 
             return(results)
            }

