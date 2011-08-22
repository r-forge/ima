lumi.methy450R <-
function(file = "../data/SampleMethFinalReport_Smiraglia.txt",
  columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = "../data/sample.txt"){
        data = read.delim(file,sep="\t",skip =8, row.names = 1)
        cat("dimension of the input methylation data",dim(data),"\n")
        #####beta matrix
        betamatrix = data[,grep(columnGrepPattern$beta,colnames(data))]
        #####detection Pvalue
        detect_p = data[,grep(columnGrepPattern$detectp,colnames(data))]
        #####Annotation
        index = grep("ILMNID",colnames(data))
        annotation = data[,index:ncol(data)]

        sname = sub(columnGrepPattern$beta,"",colnames(betamatrix))
        colnames(betamatrix) = colnames(detect_p) = sname
        group = read.table(groupfile,sep="\t")
        df = data.frame(samplename=as.character(group[,1]),groupinfo = group[,2],row.names = sname)
         setClass("exprmethy450", representation(bmatrix = "matrix",annot = "matrix",detectP = "matrix",groupinfo = "list"))
        x.methy450 = new("exprmethy450",bmatrix=as.matrix(betamatrix),annot=as.matrix(annotation),detectP=as.matrix(detect_p),groupinfo = df)
        cat("Slot names of x.methy450:",slotNames(x.methy450),"\n")
        return(x.methy450)
}

