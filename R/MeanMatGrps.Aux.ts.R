MeanMatGrps.Aux.ts <-
function(mat.list,num.time.pnts)
                {
                 I=length(mat.list)
                 results=matrix(0,num.time.pnts,num.time.pnts)	
                 for(i in 1:I)
                    results=results+mat.list[[i]]
                 results=results/I
                 results   	
                }

