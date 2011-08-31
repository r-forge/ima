testfunc <-
function(eset=eset,concov=c("ON","OFF"),testmethod = c("wilcox","limma","pooled","satterthwaite"),Padj=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none"),grouplev =grouplev,gcase = "g1",gcontrol="g2"){
        if(concov == "ON"){
              cat("RUN Robust linear regression\n")
              require(MASS)
              covariate = as.numeric(as.character(grouplev))
              testout = apply(eset,1,function(x){temp =summary(rlm(x~grouplev));pvalue = pt(abs(temp$coefficients[2,3]),4,lower.tail=FALSE)*2;return(pvalue)});
          }else{
                lev1 = which(grouplev == gcase)
                lev2 = which(grouplev == gcontrol)
                eset = eset[,c(lev1,lev2)]
                if(testmethod == "wilcox"){
                       cat("wilcox\n")
                       testout = apply(eset,1,function(x){wilcox.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))])$p.value})
                }
                if(testmethod =="limma"){
                        cat("RUN Limma\n")
                        require(limma);
                        TS = as.factor(c(rep("case",length(lev1)),rep("control",length(lev2))));
                        design = model.matrix(~0+TS);
                        rownames(design) = colnames(eset);
                        colnames(design) = c("case","control");
                        print(design);
                        fit = lmFit(eset, design);
                        cont.matrix = makeContrasts(
                                comp = case - control,
                                levels=design);
                        fit2 = contrasts.fit(fit, cont.matrix);
                        fit2 = eBayes(fit2);
                        result1 = topTable(fit2,coef=1, adjust=Padj,number=nrow(fit2));
                        testout = result1[match(rownames(eset),result1[,1]),"P.Value"];
                      }
                if(testmethod == "pooled"){
                        cat("RUN pooled t-test\n")
                        testout = apply(eset,1,function(x){t.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))],var.equal = TRUE)$p.value})
                 }
                if(testmethod == "satterthwaite"){
                        cat("RUN satterthwaite\n;")
                        testout = apply(eset,1, function(x){t.test(x[1:length(lev1)],x[(length(lev1)+1):(length(lev1)+length(lev2))])$p.value});
                }
          }
                adjustP = p.adjust(testout,method = Padj)
                if(concov=="OFF"){
                difb = apply(eset,1, function(x){mean(x[(length(lev1)+1):(length(lev1)+length(lev2))])-mean(x[1:length(lev1)])});
                out = cbind(testout,adjustP, difb,eset);
                rownames(out) = rownames(eset);
                colnames(out) = c("P-Value", "Adjust Pval","beta-Difference",colnames(eset));
              }else{
                out = cbind(testout,adjustP,eset);
                rownames(out) = rownames(eset);
                colnames(out) = c("P-Value", "Adjust Pval",colnames(eset));
              }
                return(out)
        }

