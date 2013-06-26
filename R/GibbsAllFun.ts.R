GibbsAllFun.ts <-
function(y.mu.a,y.mu.b,grp.sz,beta.mat,alfa.mat,sgm.alfa,rho,pi.i,mm,aa.pi,lmbd,df.lmbd,lmbd.alf,indexes,num.time.pnts,num.pcnts,apriori.diff.exp,nsim,burn.in,often,prob.cut.off)
              {
               if(indexes[1]!=0)
               {warning("indexes that delinatethe different time points for the sample identifiers need to begin in 0. See the manual of the function for help.")}
               if(nsim<=burn.in)
               {stop("The number of simulations is less than the burn-in period. See the package vignette for help.")}  
               nn=length(grp.sz)
               pb=tcltk::tkProgressBar(title="Gibb's progress bar", min=0,max=nsim, width=300)	
               for(i in 2:nsim)
                  {
                   results=matrix(0,nn*num.time.pnts,num.time.pnts)
                   results.SgmY=SgmYCFun.ts(nn,num.time.pnts,num.pcnts,grp.sz,y.mu.a,y.mu.b,beta.mat[,,(i-1)],alfa.mat[,,(i-1)],lmbd*df.lmbd,df.lmbd,indexes,results)
                   
                   results=matrix(0,nn,num.time.pnts)
                   beta.mat[,,i]=BetaCFun.ts(nn,num.time.pnts,num.pcnts,grp.sz,y.mu.a,y.mu.b,alfa.mat[,,(i-1)],results.SgmY,indexes,results)
                   
                   results=matrix(0,nn,num.time.pnts)
                   alfa.mat[,,i]=AlfaAllCFun.ts(nn,num.time.pnts,num.pcnts,grp.sz,y.mu.b,beta.mat[,,i],results.SgmY,sgm.alfa[,,1],rho[(i-1)],mm,indexes,results)
                   
                   pi.i[,i]=PiFun(alfa.mat[,1,i],rho[(i-1)],aa.pi,mm)
	               sgm.alfa[,,2]=SgmAlfaFun.ts(alfa.mat[,,i],lmbd.alf*num.time.pnts,num.time.pnts)
                   rho[i]=RhoFun(pi.i[,i],nn,(apriori.diff.exp/nn))
                   
                   sgm.alfa[,,1]=sgm.alfa[,,2]
                   if(i%%often==0)
                      {
                       print(i)
                       if(i>burn.in)
                         {
                          total=length(beta.mat[1,1,burn.in:i])
                          medias.pi.a=apply(alfa.mat[,1,burn.in:i],1,function(x) length(which(x!=0))/total)
                          selec=which(medias.pi.a>=prob.cut.off)
                          print(length(selec))	
                         }
                       }
                   tcltk::setTkProgressBar(pb,i,label=paste(round(i/nsim*100,0),"% done"))       
	              }
	            close(pb)  
	            resultado=list(beta.mat,alfa.mat,pi.i,rho)
                names(resultado)=c("Beta","Alfa","Pi","Rho")
                return(resultado)  
	           }

