MakeAffyOrderGrp.Fun <-
function(gene.GO.grp,affy.indexes,genes)
                    {
                     comparison=match(genes,unlist(gene.GO.grp),nomatch=0)
                     x=affy.indexes[comparison!=0]
					 x
					}

