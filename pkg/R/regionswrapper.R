regionswrapper <-
function(dataf,sumregion=c("mean","median","tbrm"),gcase = "g2",gcontrol="g1",testmethod = c("wilcox","limma","pooled","satterthwaite"),Padj="BH",concov = c("OFF","ON"),list14excel= list14excel,list14Rdata = list14Rdata){
        require("dplR")
        beta = dataf@bmatrix;
        group = dataf@groupinfo;
        index = match(colnames(beta),group$samplename)
        cat("Make sure the samples' name in group and samples' name in beta matrix are matched\n",index,"\n")
        grouplev = group$group[index]
        if(concov=="OFF"){
                lev1 = which(grouplev == gcase)
                lev2 = which(grouplev == gcontrol)
                beta = beta[,c(lev1,lev2)]
        }
         list14 = c("TSS1500Ind","TSS200Ind","UTR5Ind","EXON1Ind","GENEBODYInd","UTR3Ind", "PROMOTERInd", "GENEInd", "ISLANDInd", "NSHELFInd","NSHOREInd","SSHELFInd", "SSHOREInd","UCSCInd");
        for(i in 1:14){
                cat("calculating",list14[i],"\n");
                eval(parse(text=paste("indexlist=dataf@",list14[i],sep="")));
                eset = sumregionfun(indexlist,beta,sumregion);
                cat(list14[i],": Sum region done\nStart testing\n");
                eval(parse(text=paste(list14[i],"test=as.data.frame(testfunc(eset,concov = concov,testmethod=testmethod,Padj=Padj,lev1=lev1,lev2=lev2,grouplev=grouplev))",sep="")));
        }

        save(TSS1500Indtest,TSS200Indtest,UTR5Indtest,EXON1Indtest,GENEBODYIndtest,UTR3Indtest,PROMOTERIndtest,GENEIndtest,ISLANDIndtest,NSHELFIndtest,NSHOREIndtest,SSHELFIndtest,SSHOREIndtest,UCSCIndtest,file = "list14Rdata")
       require(WriteXLS);
       WriteXLS(paste(list14,"test",sep=""),ExcelFileName = list14excel, SheetNames = list14,row.names=TRUE)
}

