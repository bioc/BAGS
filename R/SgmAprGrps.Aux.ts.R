SgmAprGrps.Aux.ts <-
function(y.mu.a,y.mu.b,indexes,num.time.pnts,num.pcnts)
                  {
                   vec.B=matrix(0,2*num.pcnts,num.time.pnts)
                   for(j in 0:(num.time.pnts-1))
                      {
                       indices=indexes+num.pcnts*j
                       vec.B[1:num.pcnts,(j+1)]=y.mu.a[indices]
                       vec.B[(num.pcnts+1):(2*num.pcnts),(j+1)]=y.mu.b[indices]
                      }
                   lambda.0=var(vec.B)
                   lambda.0
                  }

