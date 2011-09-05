sumregionfun <-
function(indexlist,beta,sumregion){
                temp2 = matrix(NA,nrow= length(indexlist),ncol = ncol(beta))
                rownames(temp2) = names(indexlist)
                colnames(temp2) = colnames(beta)
                for(i in 1:length(indexlist)){
                       require("dplR")
                        #if(i%%1000==0){cat(i)};
                        temp = beta[indexlist[[i]],]
                        if(length(indexlist[[i]]) ==1){
                                temp2[i,] = temp;
                        }else{
                                temp2[i,] = apply(temp,2,eval(sumregion),na.rm=TRUE)
                        }
                        
                  }
                 return(temp2)
        }

