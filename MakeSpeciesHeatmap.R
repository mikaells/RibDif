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
species = args[2] # e.g subtilis
outPath = str_replace(string = ucPath,pattern = ".uc",replacement = paste("-",species,".pdf",sep=""))

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
  ucFile=rbind(ucFile[-which(ucFile$species=="sp."),],ucFile[which(ucFile$species=="sp."),])
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
#colnames(clusterMat)=unique(ucFile$V2)
colnames(clusterMat)=0:max(unique(ucFile$V2))

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
  clustMatch=which(clusterMat[,i]>0)
  #find unique species in this cluster
  unqSpecCluster=unique(rowAnnotVec[clustMatch])
  #if there are more than 1 species in the cluster
  if(length(unqSpecCluster)>1) {
    #collapse these species by '/' and add to combinations
    combinations=c(combinations,paste(sort(unqSpecCluster), collapse = "/"))
  }
}




#####
#Prepare metadata for heatmaps
#####

#build rowAnnot df as the species of each GCF
rowAnnot=data.frame(Species=rowAnnotVec)
#match rownames as GCF
rownames(rowAnnot)=unique(ucFile$GCF)

#####
#This is where it differs from standard analysis
#####
#subset for individual species

speciestIndx=rowAnnot$Species==species
speciesAlleIndx=which(colSums(clusterMat[speciestIndx,])>0)
commonSpeciesGCFIndx=which(rowSums(clusterMat[,speciesAlleIndx])>0)


clusterMat=clusterMat[commonSpeciesGCFIndx,speciesAlleIndx]



rowAnnot=rowAnnot[commonSpeciesGCFIndx,,drop=F]

#make colorfunction for species colors
colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue", "black"))
#make colors for species 
cols=colfunc(length(unique(rowAnnot$Species)))

#make a list for species/colors heatmap
annot_cols=list(Species=cols[factor(unique(rowAnnot$Species))])
names(annot_cols$Species)=unique(rowAnnot$Species)

#overwrite '.sp' as grey
annot_cols$Species[which(names(annot_cols$Species)=="sp.")]="grey"

#create cluster objects for heatmap on binary data
binClusterMat=ifelse(clusterMat>0,1,0)

colClusterMat=hclust(dist(t(binClusterMat)),"ward.D2")
rowClusterMat=hclust(dist(binClusterMat),"ward.D2")

#####
#Make heatmaps
#####

#Heatmap for cluster distribution

pdf(outPath, onefile=T)
setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=0.99, name="vp", just=c("right","top"))), action="prepend")



pheatmap(clusterMat,cluster_cols = colClusterMat,cluster_rows = rowClusterMat,  
         fontsize_col = 6,fontsize_row = 6, fontsize = 8,
         display_numbers = F,annotation_colors =annot_cols,labels_row = rowAnnot$Species,  
         number_format = "%.0f", annotation_row = rowAnnot, 
         clustering_method = "ward.D2" )

setHook("grid.newpage", NULL, "replace")
grid.text("Allele index", y=0.045,x=.4, gp=gpar(fontsize=10, fontface="bold"))

dev.off()

