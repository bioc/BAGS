GrpMeanDataFun <-
function(data.grps,data,phenotype.ids)
                {	
                 affy.ids=unlist(data.grps)                 
                 results=apply(data[affy.ids,phenotype.ids],2,mean)
                 results
                }

