IMA.methy450PP <-
function(data,na.omit = TRUE,normalization=FALSE,transfm = FALSE,samplefilterdetectP = c(FALSE,1e-5),samplefilterperc = 0.75,sitefilterdetectP = c(FALSE,0.05),sitefilterperc=0.75,locidiff = c(FALSE,0.01), Xchrom = TRUE){
	bmatrix = data@bmatrix
        detect_p = data@detectP
        annotation = data@annot;
        groupinfo = data@groupinfo
        orignalrownm = rownames(bmatrix)
        control = groupinfo$samplename[groupinfo$group=="g1"]
        experiment = groupinfo$samplename[groupinfo$group=="g2"]
        if(samplefilterdetectP){
             goodsample   = colSums(detect_p<=samplefilterdetectP)>=samplefilterperc*nrow(detect_p)
             bmatrix = bmatrix[,goodsample]
             detect_p = detect_p[,goodsample]
             cat("Total samples:",colnames(bmatrix),"\n","Kept samples:",colnames(bmatrix)[goodsample],"\n");
             groupinfo =list(samplename = as.character(groupinfo$samplename[goodsample]),group=groupinfo$group[goodsample])
              control = groupinfo$samplename[groupinfo$group=="g1"]
              experiment = groupinfo$samplename[groupinfo$group=="g2"]
        }
        if(na.omit){###Remove the sites containing missing value
              bmatrix = na.omit(bmatrix);
              temp = orignalrownm%in%rownames(bmatrix);
              detect_p = detect_p[temp,];
              annotation = annotation[temp,];
              temp = nrow(data@bmatrix)-nrow(bmatrix);
              cat(temp,"sites containing missing value and removed","\n");
         }

        if(Xchrom){
                chr = annotation[,"CHR"]
                index = which(chr == "X") #Filtering  the loci from chromsome X, total is 11232 on illumina infinium human 450 beadchip
                good_chrom = rownames(annotation)[-index]
                cat(length(index),"sites on chrX and removed\n")

        }else{good_chrom = rownames(bmatrix)}
        if(sitefilterdetectP){
           good_loci = rownames(detect_p)[rowSums(detect_p<=sitefilterdetectP)>=sitefilterperc*ncol(detect_p)]
         }else{good_loci = rownames(bmatrix)}
        if(normalization){###quantile normalization
              require(preprocessCore)
              bmatrix = normalize.quantiles(as.matrix(bmatrix))
              colnames(bmatrix) = colnames(detect_p)
              rownames(bmatrix) = rownames(detect_p)
              cat("quantile normalization Performed\n")
        }
        if(transfm){###transfer beta matrix by the arcsin square root 
                if(na.omit){bmatrix = asin(sqrt(bmatrix));cat("transfer beta matrix by the arcsin square root\n")}
                else{cat("\tMissing value exist in the orignial data,\nPlease remove the missing value before transformation,use na.omit = TRUE\n");
                stop;}
       }
        if(locidiff){
                con_mean = apply(bmatrix[,control],1,mean)
                trt_mean = apply(bmatrix[,experiment],1,mean)
                good_diff = rownames(bmatrix)[abs(1-trt_mean/con_mean)>=locidiff]

        }else{
                good_diff = rownames(bmatrix)
        }
        all_good = intersect(intersect(good_chrom,good_loci),good_diff);
        cat(length(all_good), "were retained from the total sites",length(orignalrownm),"\n")
        bmatrix = bmatrix[all_good,]
        annotation = annotation[all_good,]
        detect_p = detect_p[all_good,]
        cat(".......\n split the annotation file to 11 annotated region list\n .......\n" )
##############################################################        
#########split the annotation file to 11 annoted region list
###Make the new annotation file by the UCSC_REFGENE_NAME
        annot = annotation
        name = "UCSC_REFGENE_NAME"
        cpGsite = as.character(annot[,1])
        genelist = strsplit(as.character(annot[,name]),";");
        genelist[which(genelist == "character(0)")]="NA"
        name = "UCSC_REFGENE_GROUP"
        refgene = strsplit(as.character(annot[,name]),";")
        refgene[which(refgene == "character(0)")] = "NA"
        listlength = lapply(refgene,length)
        listlength[listlength ==0] = 1
        col1 = rep(cpGsite, listlength)##cpGsite
        col2 = unlist(genelist)##GeneSymbol
        col3 = unlist(refgene)##REFGENE_GROUP
        col4 = rep(as.character(annotation[,"RELATION_TO_UCSC_CPG_ISLAND"]),listlength) ###RELATION_TO_UCSC_CPG_ISLAND
        col5 = rep(as.character(annotation[,"UCSC_CPG_ISLANDS_NAME"]),listlength)
###UCSC_CPG_ISLANDS_NAME
        annotation2 = cbind(col1,col2,col3,col4,col5);
        colnames(annotation2) = c("cpGsite","GeneSymbol","REFGENE_GROUP","RELATION_TO_UCSC_CPG_ISLAND","UCSC_CPG_ISLANDS_NAME");
####split the annotation to by the REF Gene Group:list 1-8
        splitToRegionlist = function(grepname = c("TSS1500","TSS200","5'UTR","1stExon","Gene Body","3'UTR")){
                index = col3 == grepname
                col1sub = col1[index]
                col2sub = col2[index]
                temp = split(col1sub,col2sub)
                returnInd = lapply(temp,unique)###may have some duplicates here
                return(returnInd)
         }
###REF Gene Group
        TSS1500Ind = splitToRegionlist(grepname="TSS1500")
        TSS200Ind = splitToRegionlist(grepname = "TSS200")
        UTR5Ind = splitToRegionlist(grepname = "5'UTR")
        EXON1Ind = splitToRegionlist(grepname = "1stExon")
        GENEBODYInd = splitToRegionlist(grepname = "Body")
        UTR3Ind = splitToRegionlist(grepname = "3'UTR")
        cat("TSS1500 region:",length(TSS1500Ind),"\nTSS200 region:",length(TSS200Ind),"\n5'UTR region:",length(UTR5Ind),"\n1st Exon region:",length(EXON1Ind),"\nGene body region:",length(GENEBODYInd),"\n3'UTR region:",length(UTR3Ind),"\n")

####Now Split by the RELATION_TO_UCSC_CPG_ISLAND
        splitToRegionlist2 = function(grepname = c("Island","N_Shore","S_Shore","N_Shelf","S_Shelf")){
                index = col4 == grepname
                col1sub = col1[index]
                col5sub = col5[index]
                temp = split(col1sub,col5sub)
                returnInd = lapply(temp,unique)###may have some duplicates here
                return(returnInd)
        }
        ISLANDInd = splitToRegionlist2(grepname="Island")
        NSHOREInd = splitToRegionlist2(grepname="N_Shore")
        SSHOREInd = splitToRegionlist2(grepname="S_Shore")
        NSHELFInd = splitToRegionlist2(grepname="N_Shelf")
        SSHELFInd = splitToRegionlist2(grepname="S_Shelf")
        cat("Island region:",length(ISLANDInd),"\nN_Shore region", length(NSHOREInd),"\nS_Shore region", length(SSHOREInd),"\nN_Shelf region",length(NSHELFInd),"\nS_Shelf region",length(SSHELFInd),"\n")
        setClass("methy450batch", representation(bmatrix = "matrix",annot = "matrix",detectP = "matrix",groupinfo = "list",TSS1500Ind = "list", TSS200Ind = "list",UTR5Ind = "list", EXON1Ind = "list", GENEBODYInd = "list",UTR3Ind = "list",ISLANDInd="list",NSHOREInd = "list",SSHOREInd="list",NSHELFInd="list",SSHELFInd ="list"))
        x.methy450 = new("methy450batch",bmatrix=as.matrix(bmatrix),annot=as.matrix(annotation),detectP=as.matrix(detect_p),groupinfo = groupinfo,TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind,UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind,GENEBODYInd = GENEBODYInd, UTR3Ind = UTR3Ind,ISLANDInd=ISLANDInd,NSHOREInd = NSHOREInd,SSHOREInd=SSHOREInd,NSHELFInd=NSHELFInd,SSHELFInd =SSHELFInd)
        cat("Slot names of methy450batch:",slotNames(x.methy450),"\n")
        return(x.methy450)
}

