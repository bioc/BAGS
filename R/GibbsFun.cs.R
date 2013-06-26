GibbsFun.cs <-
function(y.mu.a,y.mu.b,grp.sz,beta.mat,alfa.mat,sgm.y.a,sgm.y.b,sgm.alfa,rho,pi.i,mm,aa,bb,aa.pi,apriori.diff.exp,nsim,burn.in,often,prob.cut.off)
             {
              if(aa<=0 | bb<=0)
              {stop("Parameters of the prior inverse gamma distribution for the variance parameter assumed on the data are not valid. They need to be both positive. See reference for more information")}	
              if(nsim<=burn.in)
               {stop("The number of simulations is less than the burn-in period. See the package vignette for help.")}
              if(mm<=0 | aa.pi<=0)
              {stop("Parameters of the prior beta distribution for the non-zero part of the parameter Pi are not valid. They need to be both positive. See reference for more information")}
              	
              pb=tcltk::tkProgressBar(title="Gibb's progress bar", min=0,max=nsim, width=300)	
              for(i in 2:nsim)
   				{
       	         sgm.y.a[,i]=SgmYFun.cs(grp.sz,y.mu.a,beta.mat[,(i-1)],bb,aa)
                 sgm.y.b[,i]=SgmYFun.cs(grp.sz,y.mu.b,beta.mat[,(i-1)]+alfa.mat[,(i-1)],bb,aa)
   	
                 beta.mat[,i]=BetaFun.cs(y.mu.a,y.mu.b,grp.sz,alfa.mat[,(i-1)],sgm.y.a[,i],sgm.y.b[,i])	
   	             alfa.mat[,i]=AlfaFun.cs(y.mu.b,grp.sz,beta.mat[,i],sgm.y.b[,i],sgm.alfa[(i-1)],rho[(i-1)],mm)
    
                 pi.i[,i]=PiFun(alfa.mat[,i],rho[(i-1)],aa.pi,mm)
                 sgm.alfa[i]=SgmAlfaFun.cs(alfa.mat[,i],aa,bb)
                 rho[i]=RhoFun(pi.i[,i],length(grp.sz),(apriori.diff.exp/length(grp.sz)))
    
                 if(i%%often==0)
                   {
                    print(i)
                    if(i>burn.in)
                      {
                       total=length(beta.mat[1,burn.in:i])
                       medias.pi.a=apply(alfa.mat[,burn.in:i],1,function(x) length(which(x!=0))/total)
                       selec=which(medias.pi.a>prob.cut.off)
                       print(length(selec))	
                      }
                   }
                 tcltk::setTkProgressBar(pb,i,label=paste(round(i/nsim*100,0),"% done"))  
                }
              close(pb)
              results=list(beta.mat,sgm.y.a,alfa.mat,sgm.y.b,pi.i,sgm.alfa,rho)  
              names(results)=c("Beta","Sgm.Y.A","Alfa","Sgm.Y.B","Pi","Sgm.alfa","Rho")
              return(results)
             }

