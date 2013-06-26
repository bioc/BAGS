GeneGrps2AffyGrpsFun<-
function(GO.group,level,gene.u.name,affy.max.indexes)
                        {
                         if(level<3 | level>10)
                         {stop("The depth of the hierarchical has to be between 3 and 10. See help for more information")}
                         if(is.element(GO.group,c("CC","MF","BP"))==FALSE)
                         {stop("The choice of ontology has to be of the form: CC (Cellular component), MF (Molecular function) or BP (Biological processes)")}
                         
                         data(GGAdata,package="GeneGroupAnalysis") 	
                         process.tree=switch(GO.group,CC=GGAdata$CClevels,MF=GGAdata$MFlevels,BP=GGAdata$BPlevels)	
              	         gene.grouping=process.tree[[level]]
              	         out.aux=lapply(1:length(gene.grouping),function(x){length(gene.grouping[[x]])})
                         out=which(unlist(out.aux)==0)
                         gene.grouping.b=gene.grouping[-out]
                         GO.process.b=names(gene.grouping.b)
                         
                         affy.grouping=lapply(1:length(gene.grouping.b),function(x) MakeAffyOrderGrp.Fun(gene.grouping.b[x],affy.max.indexes,gene.u.name))
                         out.aux=lapply(1:length(affy.grouping),function(x){length(affy.grouping[[x]])})
                         out=which(unlist(out.aux)==0)
                       
                         results=list(affy.grouping[-out],gene.grouping.b[-out],GO.process.b[-out])
                         names(results)=c("index.GO.grps","gene.GO.grps","GO.grps")
                         return(results)
                        }
