#####
#Setup libs
#####
suppressMessages(library(ape, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(stringr, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(adephylo, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(pheatmap, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(grid, quietly = T, warn.conflicts = FALSE))


#####
#handle input arguments
#####
args=commandArgs(trailingOnly = TRUE);

#get command args
ucPath  = args[1] #vsearch clusters (.uc-file)
outPath = args[2] #path for pdf etc


#####
#Read and clean-up .uc-file
#####

#read cluster file
ucFile=read.table(ucPath)

#removing cluster summaries
ucFile=ucFile[-grep("C",ucFile$V1),]

#get GCF and species from names
ucFile$GCF=apply(str_split_fixed(ucFile$V9,"_", 8)[,1:2] , 1 , paste , collapse = "_" )
ucFile$species=str_split_fixed(ucFile$V9,"_", 8)[,6]

#sort UC by GCF to keep genomes together and reproducible
ucFile=ucFile[order(ucFile$species,ucFile$GCF),]

#putting all unclassified species at the bottom
if(length(which(ucFile$species=="sp."))) {
  ucFile=rbind(ucFile[-c(ucFile$species=="sp."),],ucFile[-c(ucFile$species=="sp."),])
}

#####
#Prepare matrix for cluster abundances
#The idea is to have a unique cluster in each column
#and then have genomes (GCFs) represented in rows
#####

#create empty matrix for storing cluster abundances
#rows are GCFs and cols are clusters
clusterMat=data.frame(matrix(0,nrow=length(unique(ucFile$GCF)),ncol = 1+length(unique(ucFile$V2))))
rownames(clusterMat)=unique(ucFile$GCF)
colnames(clusterMat)=unique(ucFile$V2)

#making an empty vector for storing species corresponding to GCFs
rowAnnotVec=rep("",length(unique(ucFile$GCF)))

#make a counter for rowIndx in matrix
rowCount=1

#loop through each GCF 
for(j in (unique(ucFile$GCF))){
  #find indeces for all 16s genes for a given GCF
  GCFindx=which(ucFile$GCF %in% j)
  #subset to only 16s genes for a given GCF
  ucClus=ucFile[GCFindx,]
  #fill in species of the GCF
  rowAnnotVec[rowCount]=unique(ucClus$species)
  #loop through all 16s genes of the GCF
  for(k in 1:NROW(ucClus)) {
    #in the row for the GCF (constant in this loop)
    #and the column corresponding to the .uc-cluster of this 16S gene
    #increment by one
    #+1 because .uc-clusters are 0-indexed
    clusterMat[rowCount,ucClus$V2[k]+1] = clusterMat[rowCount,ucClus$V2[k]+1]+1
  }
  #increment rowcount for next GCF
  rowCount=rowCount+1
}

#####
#Find all instances of species overlap
#####

#make empty combination vector
combinations=c()


#loop throug all columns of clusterMat, e.g. all .uc-clusters
for(i in 1:NCOL(clusterMat)) {
  #find all row indeces (GCFs) which have a least 1 member in cluster
  clustMatch=which(clusterMat[,i]>1)
  #find unique species in this cluster
  unqSpecCluster=unique(rowAnnotVec[clustMatch])
  #if there are more than 1 species in the cluster
  if(length(unqSpecCluster)>1) {
    #collapse these species by '/' and add to combinations
    combinations=c(combinations,paste(sort(unqSpecCluster), collapse = "/"))
  }
}

#####
#Make a matrix of all GCF overlaps
#####

#prepare empty matrix, symmetrical with rows/cols for GCF
pairwiseMatch=matrix(data=0,length(unique(ucFile$GCF)),length(unique(ucFile$GCF)))
#give row and colnames as GCF
colnames(pairwiseMatch)= rownames(pairwiseMatch)=unique(ucFile$GCF)

j=unique(ucFile$GCF)[1]
#make rowCount for row bookkeeping
rowCount=1

#for each unique GCF
for(j in (unique(ucFile$GCF))){
  #visual counter
  #if(rowCount %in% seq(1, 1700,by = 50)) print(rowCount)
  #find indeces for GCF
  GCFindx=which(ucFile$GCF %in% j)
  #subset to GCF
  ucClus=ucFile[GCFindx,]
  #find species for GCF
  clusSpec=unique(ucClus$species)
  #find unique clusters in GCF
  clusClus=unique(ucClus$V2)
  #find the GCFs that have genes in GCF of the loop
  clusMatchGCF=unique(ucFile$GCF[which(ucFile$V2 %in% clusClus)])
  
  #for each of these matching GCFs
  for(k in clusMatchGCF) {
    #get the index of the matching GCF (column)
    colIndx=which(colnames(pairwiseMatch) %in% k)
    #and add a 1 to the row of original GCF and col of new GCF
    pairwiseMatch[rowCount,colIndx]=1
  }
  #update row/GCF
  rowCount=rowCount+1
}


#####
#Prepare metadata for heatmaps
#####

#build rowAnnot df as the species of each GCF
rowAnnot=data.frame(Species=rowAnnotVec)
#match rownames as GCF
rownames(rowAnnot)=unique(ucFile$GCF)

#make colorfunction for species colors
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue", "black"))
#make colors for species 
cols=colfunc(length(unique(rowAnnot$Species)))

#make a list for species/colors heatmap
annot_cols=list(Species=cols[factor(unique(rowAnnot$Species))])
names(annot_cols$Species)=unique(rowAnnot$Species)

#overwrite '.sp' as grey
annot_cols$Species[which(names(annot_cols$Species)=="sp.")]="grey"

#####
#Make heatmaps
#####

#Heatmap for cluster distribution

pdf(outPath, onefile=T)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=0.99, name="vp", just=c("right","top"))), action="prepend")

pheatmap(clusterMat, cluster_rows = T, fontsize_col = 6,fontsize_row = 8, fontsize = 8,
         cluster_cols = T, display_numbers = F,annotation_colors =annot_cols,labels_row = rowAnnot$Species,  
         number_format = "%.0f", annotation_row = rowAnnot, clustering_method = "ward.D2" )

setHook("grid.newpage", NULL, "replace")
grid.text("Allele index", y=0.045,x=.4, gp=gpar(fontsize=10, fontface="bold"))

#Heatmap for pairwiseMatch 

pheatmap(pairwiseMatch,  annotation_colors =annot_cols, fontsize = 8, 
         annotation_row = rowAnnot, labels_row = rowAnnotVec, clustering_method = "ward.D2")
		   
garbage <- dev.off()

#write out table
write.table(x = pairwiseMatch,file = gsub("pdf","csv",outPath),
            row.names = rowAnnotVec,sep = ";")

#####
#Tabulate data and print 
#####

#Total genomes
TotGCF=NROW(clusterMat)
#n genomes with multiple alleles
multiAllele=length(which(rowSums(ifelse(clusterMat>0,1,0))>1))

if(length(unique(combinations))>0){
  max_nOverlaps=max(str_count(combinations,pattern = "/")+1)
  #counting species that overlap by counting '/'
  hasOverlap=unique(as.character(str_split_fixed(combinations,pattern = "/",n = max_nOverlaps)))
} else {
  max_nOverlaps=0
  #counting species that overlap by counting '/'
  hasOverlap=unique(as.character(str_split_fixed(combinations,pattern = "/",n = max_nOverlaps)))
}

#removing 'sp'
hasOverlap=hasOverlap[-grep("^$|sp\\.",hasOverlap)]

#counting overlapping species
overlapSpecN=length(hasOverlap)
#counting total species
totSpec=length(unique(rowAnnot$Species[!rowAnnot$Species=="sp."]))

#total named species
nNamedGCF=length((rowAnnot$Species[!rowAnnot$Species=="sp."]))
#total non-named species
nNonnamedGCF=length((rowAnnot$Species[rowAnnot$Species=="sp."]))

#making summary string
summaryString=paste("Summary:\n\n",
                    "Genomes: ",TotGCF,"\n",paste("\t","Named: ",nNamedGCF,"\n\tNon-named: ",nNonnamedGCF,"\n\n", sep=""),
                    "Named species: ",totSpec,"\n\n",
                    multiAllele, " of ", TotGCF, " (",signif(100*multiAllele/TotGCF,4),"%) genomes have multiple alleles.\n\n",
                    overlapSpecN," of ", totSpec," (",signif(100*overlapSpecN/totSpec,4),"%) species overlap.\n\n",
                    ifelse(length(unique(combinations))>0,paste("The following species overlap:\n\t",paste(unique(combinations),"\n", collapse = "\t",sep="")),
                           paste("No species overlap.\n"))
                    ,sep="")

#print to terminal
cat(summaryString)

#add to log file
sink(file =  gsub("heatmap.pdf","overlap-summary.txt",outPath))
cat(summaryString)
sink()

