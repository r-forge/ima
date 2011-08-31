IMA.methy450R <-
function(file = MethtFileName, columnGrepPattern=list(beta=".AVG_Beta",detectp=".Detection.Pval"),groupfile = PhenoFileName){
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
        ####match the phenotype data with the methylation data
        group = read.table(groupfile,sep="\t",header = TRUE)
        index = match(group[,1],colnames(betamatrix))
        betamatrix=betamatrix[,index]
        detect_p = detect_p[,index]
        df = data.frame(samplename=as.character(group[,1]),groupinfo = group[,2],row.names = sname)
         setClass("exprmethy450", representation(bmatrix = "matrix",annot = "matrix",detectP = "matrix",groupinfo = "list"))
        x.methy450 = new("exprmethy450",bmatrix=as.matrix(betamatrix),annot=as.matrix(annotation),detectP=as.matrix(detect_p),groupinfo = df)
        cat("Slot names of x.methy450:",slotNames(x.methy450),"\n")
        require(bioDist)
        eset = na.omit(betamatrix)
        samples = paste(group[,1],group[,2],sep="_")
        hc1<-hclust(cor.dist(t(eset)), method="average")
        pdf("./QC.pdf")
        #postscript("../../pg_0001.ps")
        plot(hc1,samples, xlab="Sample", main="Clusting All Genes (Pearson correlation)", lwd=2, font.axis=2, font.lab=2)
        #dev.off()
        #postscript("../../pg_0002.ps")
        boxplot(betamatrix,ylab = "beta Value")
        #dev.off()
        #postscript("../../pg_0003.ps")
        avgPval = apply(detect_p,2,function(x){sum(x>=1e-5)*100/length(x)})
        barplot(avgPval, ylab = "% of detect pvalue >1e-5")
        dev.off()
        return(x.methy450)
}

