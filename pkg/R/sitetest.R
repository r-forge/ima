sitetest <-
function(dataf,gcase = "g2",gcontrol="g1",testmethod = c("wilcox","limma","pooled","satterthwaite"),Padj="BH",concov = "OFF",rawpcut = NULL,adjustpcut =NULL,betadiffcut = NULL){
        beta = dataf@bmatrix;
        group = dataf@groupinfo;
        index = match(colnames(beta),group$samplename)
        cat("Make sure the samples' name in group and samples' name in beta matrix are matched\n",index,"\n")
        grouplev = group$group[index]
        lev1 = which(grouplev == gcase)
        lev2 = which(grouplev == gcontrol)
        eset = beta[,c(lev1,lev2)]###case in the first and control followed
        if(concov == "ON"){
              require(MASS)
              covariate = as.numeric(as.character(grouplev))
              testout = apply(eset,1,function(x){temp =summary(rlm(x~grouplev));pvalue = pt(abs(temp$coefficients[2,3]),4,lower.tail=FALSE)*2;return(pvalue)})

        }
        if(testmethod == "wilcox"&concov == "OFF"){
                testout = apply(eset,1,function(x){wilcox.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))])$p.value})
        }
        if(testmethod == "limma"&concov == "OFF"){
                require(limma)
                TS = as.factor(c(rep("case",length(lev1)),rep("control",length(lev2))))
                design = model.matrix(~0+TS)
                rownames(design) = colnames(eset)
                colnames(design) = c("case","control")
                fit = lmFit(eset, design)
                cont.matrix = makeContrasts(
                        comp = case - control,
                        levels=design)
                fit2 = contrasts.fit(fit, cont.matrix)
                fit2 = eBayes(fit2)
                result1 = topTable(fit2,coef=1, adjust=Padj,number=nrow(fit2))
                testout = result1[match(rownames(beta),result1[,1]),"P.Value"]
        }
        if(testmethod == "pooled"&concov == "OFF"){
                testout = apply(eset,1,function(x){t.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))],var.equal = TRUE)$p.value})
        }
        if(testmethod == "satterthwaite"&concov == "OFF"){
                testout = apply(eset,1, function(x){t.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))])$p.value})

        }
        adjustP = p.adjust(testout,method = Padj)
        difb = apply(eset,1, function(x){mean(x[1:length(lev1)])-mean(x[(length(lev1)+1):(length(lev1)+length(lev2))])})
        out = cbind(testout,adjustP, difb)
        rownames(out) = rownames(eset);
        colnames(out) = c("P-Value", "Adjust Pval","beta-Difference")
        out = outputDESfunc(out = out,rawpcut = rawpcut,adjustpcut =adjustpcut,betadiffcut = betadiffcut)
        return(out)
}

