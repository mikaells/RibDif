suppressMessages(library(ape, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(stringr, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(adephylo, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(pheatmap, quietly = T, warn.conflicts = FALSE))

args=commandArgs(trailingOnly = TRUE);

#get command args
ucPath  = args[1]
outPath = args[2]


#read cluster file
ucFile=read.table(ucPath)

#removing cluster summaries
ucFile=ucFile[-grep("C",ucFile$V1),]

#get GCF and species from names
ucFile$GCF=apply(str_split_fixed(ucFile$V9,"_", 8)[,1:2] , 1 , paste , collapse = "_" )
ucFile$species=str_split_fixed(ucFile$V9,"_", 8)[,6]

#create empty matrix for storing cluster abundances
#rows are GCFs and cols are clusters
clusterMat=data.frame(matrix(0,nrow=length(unique(ucFile$GCF)),ncol = length(unique(ucFile$V2))))
rownames(clusterMat)=unique(ucFile$GCF)
colnames(clusterMat)=unique(ucFile$V2)

#making an empty vector for storing species corresponding to gcfs
rowAnnotVec=rep("",length(unique(ucFile$GCF)))

rowCount=1
#Turning cluster file into a matrix 
for(j in (unique(ucFile$GCF))){
  
  GCFindx=grep(j,ucFile$GCF)
  ucClus=ucFile[GCFindx,]
  rowAnnotVec[rowCount]=unique(ucClus$species)
  for(k in 1:NROW(ucClus)) {
    clusterMat[rowCount,ucClus$V2[k]+1] = clusterMat[rowCount,ucClus$V2[k]+1]+1
  }
  
  rowCount=rowCount+1
}


i=3
combinations=c()
for(i in 1:NCOL(clusterMat)) {
  
  clustMatch=which(clusterMat[,i]>1)
  if(length(unique(rowAnnotVec[clustMatch]))>1) {
    combinations=c(combinations,paste(sort(unique(rowAnnotVec[clustMatch])), collapse = "/"))
  }
}

pairwiseMatch=matrix(data=0,length(unique(ucFile$GCF)),length(unique(ucFile$GCF)))
#rownames(pairwiseMatch)=rowAnnotVec
#colnames(pairwiseMatch)= unique(ucFile$GCF)
colnames(pairwiseMatch)= rownames(pairwiseMatch)=unique(ucFile$GCF)

j=unique(ucFile$GCF)[1]
rowCount=1
for(j in (unique(ucFile$GCF))){
  
  GCFindx=grep(j,ucFile$GCF)
  ucClus=ucFile[GCFindx,]
  clusSpec=unique(ucClus$species)
  
  clusClus=unique(ucClus$V2)
  clusMatchGCF=unique(ucFile$GCF[which(ucFile$V2 %in% clusClus)])
  
  
  for(k in clusMatchGCF) {
    
    colIndx=grep(k, colnames(pairwiseMatch))
    pairwiseMatch[rowCount,colIndx]=1
    
  }
  rowCount=rowCount+1
}


rowAnnot=data.frame(Species=rowAnnotVec)
rownames(rowAnnot)=unique(ucFile$GCF)
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue", "black"))

cols=colfunc(length(unique(rowAnnot$Species)))

annot_cols=list(Species=cols[factor(unique(rowAnnot$Species))])
names(annot_cols$Species)=unique(rowAnnot$Species)

pdf(outPath, onefile=T)

pheatmap(clusterMat, cluster_rows = T, 
                   cluster_cols = F, display_numbers = F,annotation_colors =annot_cols,  
                   number_format = "%.0f", annotation_row = rowAnnot )

pheatmap(pairwiseMatch,  annotation_colors =annot_cols,  
          annotation_row = rowAnnot, labels_row = rowAnnotVec)
				   
garbage <- dev.off()

write.table(x = pairwiseMatch,file = gsub("pdf","csv",outPath),col.names = c("",colnames(pairwiseMatch)),
            row.names = rowAnnotVec,sep = ";")

TotGCF=NROW(clusterMat)
multiAllele=length(which(rowSums(ifelse(clusterMat>0,1,0))>1))

{
cat(paste("Summary:\n"))
cat(paste(multiAllele, " of ", TotGCF, " (",signif(100*multiAllele/TotGCF,4),"%) genomes have multiple alleles.\n\n", sep=""))
if(length(unique(combinations))>0) {
  cat(paste("The following species overlap:\n"))
  cat(paste("\t",unique(combinations),"\n"))
} else {
  cat(paste("No species overlap.\n"))
}
}
