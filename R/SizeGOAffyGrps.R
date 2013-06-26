SizeGOAffyGrps <-
function(affy.GO.grps,size)
             {
              grps.sz=lapply(affy.GO.grps,function(x) length(x))
              wrkng.grps=which(unlist(grps.sz)>=size)
              results=list(wrkng.grps,unlist(grps.sz)[wrkng.grps])
              names(results)=c("groups","group.size")
              results
             }

