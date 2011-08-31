outputDESfunc <-
function(out=out,outputDES = TRUE,rawpcut = 0.05,adjustpcut =0.05,betadiffcut = 0.14){
       if(outputDES==TRUE){
          if(is.null(rawpcut)&is.null(adjustpcut)& is.null(betadiffcut)){
            cat("Please choose the cut off of the raw pvalue, adjust pvalue,or beta difference")
          }else{
                 if(!is.null(rawpcut)){
                        rawpcutout = out[,1]<=rawpcut
                  }else{
                        rawpcutout = out[,1]<=1
                  }
                  if(!is.null(adjustpcut)){
                        adjustpcutout = out[,2]<= adjustpcut
                  }else{
                       adjustpcutout = out[,2]<=1
                  }
                  if(!is.null(betadiffcut)){
                       betadiffcutout = abs(out[,3])>=betadiffcut;
                  }else{
                       betadiffcutout = abs(out[,3])>=0;
                  }
          }
        out = out[rawpcutout&adjustpcutout&betadiffcutout,]
        }else{out = out}
  }

