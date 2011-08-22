lumi.methy450PP <-
function(data,na.omit = TRUE,normalization=FALSE,transfm = FALSE,filtersample = TRUE,filterdetectP = TRUE, detectPcut = 0.75,locidiff = FALSE,locidiffcut=0.01, Xchrom = TRUE){
        bmatrix = data@bmatrix
        detect_p = data@detectP
        annotation = data@annot;
        groupinfo = data@groupinfo
        orignalrownm = rownames(bmatrix)
        control = groupinfo$samplename[groupinfo$group=="g1"]
        experiment = groupinfo$samplename[groupinfo$group=="g2"]
        if(filtersample){
             goodsample   = colSums(detect_p<= 1e-5)>detectPcut* nrow(detect_p)
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
              cat(temp,"sites containing missing value and be removed","\n");
         }

        if(Xchrom){
                chr = annotation[,"CHR"]
                index = which(chr == "X") #Filtering  the loci from chromsome X, total is 11232 on illumina infinium human 450 beadchip
                good_chrom = rownames(annotation)[-index]
                cat(length(index),"sites on chrX \n")

        }else{good_chrom = rownames(bmatrix)}
        if(filterdetectP){
                #test1 = rownames(detect_p)[rowSums(detect_p[,control]<0.01)>detectPcut*length(control)]
                #test2 = rownames(detect_p)[rowSums(detect_p[,experiment]<0.01)>detectPcut*length(experiment)]
                rowmedians = apply(detect_p,1,median)
                good_loci = rownames(detect_p)[rowmedians<=0.05]
        }else{good_loci = rownames(bmatrix)}
        if(normalization){###quantile normalization
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
                good_diff = rownames(bmatrix)[abs(1-trt_mean/con_mean)>=locidiffcut]

        }else{
                good_diff = rownames(bmatrix)
        }
        all_good = intersect(intersect(good_chrom,good_loci),good_diff);
        cat("Kept loci", length(all_good),"from orignial ",length(orignalrownm),"\n")
        bmatrix = bmatrix[all_good,]
        annotation = annotation[all_good,]
        detect_p = detect_p[all_good,]
        cat(".......\n split the annotation file to 14 annotated region list\n .......\n" )
##############################################################        
#########split the annotation file to 14 annoted region list
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
        splitToRegionlist = function(grepname = c("TSS1500","TSS200","5'UTR","1stExon","Gene Body","3'UTR","Promoter")){
                if(grepname == "Promoter"){
                        index = col3 == "TSS1500"|col3 == "TSS200"|col3 == "5'UTR"|col3 == "1stExon"
               }else if(grepname == "Gene"){
                        index = 1:length(col1)
               }else{
                        index = col3 == grepname
               }
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
        PROMOTERInd = splitToRegionlist(grepname = "Promoter")
        GENEInd = splitToRegionlist(grepname = "Gene")
        cat("TSS1500 region:",length(TSS1500Ind),"\nTSS200 region:",length(TSS200Ind),"\n5'UTR region:",length(UTR5Ind),"\n1st Exon region:",length(EXON1Ind),"\n gene body region:",length(GENEBODYInd),"\n3'UTR region:",length(UTR3Ind),"\nPROMOTER region:",length(PROMOTERInd),"\n")
#TSS1500Ind 20406 TSS200Ind 17731 UTR5Ind #14148 EXON1Ind#15588 UTR3Ind #13074 PROMOTERInd #20972
####Now Split by the RELATION_TO_UCSC_CPG_ISLAND
        splitToRegionlist2 = function(grepname = c("Island","N_Shore","S_Shore","N_Shelf","S_Shelf","UCSC")){
                if(grepname == "UCSC"){
                        index = col4 == "Island"|col4 == "N_Shelf"|col4 == "N_Shore"|col4 == "S_Shelf"|col4 == "S_Shore"

                }else{
                        index = col4 == grepname

                }
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
        UCSCInd = splitToRegionlist2(grepname="UCSC")
        cat("Island region:",length(ISLANDInd),"\nN_Shore region", length(NSHOREInd),"\nS_Shore region", length(SSHOREInd),"\n N_Shelf region",length(NSHELFInd),"\n S_Shelf region",length(SSHELFInd),"\n UCSC island region", length(UCSCInd),"\n")
 # lenght(ISLANDInd),26662  length(NSHOREInd)#[1] 24991  length(SSHOREInd)#[1] 22444  length(NSHELFInd)#[1] 18417 length(SSHELFInd)#[1] 17337 length(UCSCInd)#[1] 27176

        setClass("methy450batch", representation(bmatrix = "matrix",annot = "matrix",detectP = "matrix",groupinfo = "list",TSS1500Ind = "list", TSS200Ind = "list",UTR5Ind = "list", EXON1Ind = "list", GENEBODYInd = "list",UTR3Ind = "list",PROMOTERInd = "list", GENEInd = "list",ISLANDInd="list",NSHOREInd = "list",SSHOREInd="list",NSHELFInd="list",SSHELFInd ="list",UCSCInd = "list"))
        x.methy450 = new("methy450batch",bmatrix=as.matrix(bmatrix),annot=as.matrix(annotation),detectP=as.matrix(detect_p),groupinfo = groupinfo,TSS1500Ind = TSS1500Ind, TSS200Ind = TSS200Ind,UTR5Ind = UTR5Ind, EXON1Ind = EXON1Ind,GENEBODYInd = GENEBODYInd, UTR3Ind = UTR3Ind,PROMOTERInd = PROMOTERInd, GENEInd = GENEInd,ISLANDInd=ISLANDInd,NSHOREInd = NSHOREInd,SSHOREInd=SSHOREInd,NSHELFInd=NSHELFInd,SSHELFInd =SSHELFInd,UCSCInd = UCSCInd)
        cat("Slot names of methy450batch:",slotNames(x.methy450),"\n")
        return(x.methy450)
}

