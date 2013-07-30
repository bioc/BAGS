DataGeneSets<-
function(output.ReadGMT,data.gene.symbols,size)
		{
			#=="input": output of the ReadGMT function & vector of gene symbols from the data & minimum grp size considered
			#=="output":
			gene.nams.unclstr=unlist(lapply(1:length(data.gene.symbols),function(x){
	                      gene.clstr=strsplit(gsub("[^[:alnum:] ]"," ",data.gene.symbols[x])," ")[[1]]
                          y=nchar(gene.clstr);z=gene.clstr[which(y!=0)];return(z)}))
		    gene.ids.unclstr=unlist(lapply(1:length(data.gene.symbols),function(x){
	                      gene.clstr=strsplit(gsub("[^[:alnum:] ]"," ",data.gene.symbols[x])," ")[[1]]
                          y=nchar(gene.clstr);z=rep(x,length(which(y!=0)));return(z)}))
			grps.ids=lapply(output.ReadGMT,function(x){
				                                   aux=match(x,gene.nams.unclstr,nomatch=0)
				                                   aux=aux[which(aux!=0)];aux=gene.ids.unclstr[aux]
				                                   return(aux)})
			grps.nms=lapply(output.ReadGMT,function(x){
				                                   aux=match(x,gene.nams.unclstr,nomatch=0)
				                                   aux=x[which(aux!=0)]
				                                   return(aux)})	                                   
			wrkng.szs=unlist(lapply(grps.ids,length))
			wrkng.szs.f=wrkng.szs[which(wrkng.szs>=size)]
			wrkng.ids.f=grps.ids[which(wrkng.szs>=size)]
			wrkng.nms.f=grps.nms[which(wrkng.szs>=size)]
			
			resultado=list(wrkng.ids.f,wrkng.nms.f,wrkng.szs.f)
			names(resultado)=c("DataGeneSetsIds","DataGeneSetsNms","Size")
			return(resultado)
		}