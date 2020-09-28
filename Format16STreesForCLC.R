suppressMessages(library(ape, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(stringr, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(adephylo, quietly = T, warn.conflicts = FALSE))

args=commandArgs(trailingOnly = TRUE);

treePath = args[1]
outPath  = args[2]

#read tree
#tree=read.tree("c:/Users/milst/Desktop/Pseudoalteromonas.tree")
tree=read.tree(treePath)


#get tips for getting genes later
tips=tree$tip.label

#fixing tip names
tips1=str_split_fixed(gsub(":","",gsub("_complete_.*","",gsub("_chromosome_.*","",tips))),pattern = "_", 12)
tips1[,5]=gsub("^(.).*","\\1.",tips1[,5])

Name    = tips
GCF     = apply( tips1[ , 1:2 ] , 1 , paste , collapse = "_" )
NZ      = apply( tips1[ , 3:4 ] , 1 , paste , collapse = "_" )
Strain  = trimws(apply( tips1[ , 5:12 ] , 1 , paste , collapse = " " ), which = "right")
Species = apply( tips1[ , 5:6 ] , 1 , paste , collapse = " " )
nGene  = gsub(".*_([0-9]{1,2}$)","\\1",tips)

#Making the data frame
df=data.frame(Name, GCF, NZ, Strain, Species, nGene)

#copying tree
tree2=tree
tree2$tip.label=1:length(tree$tip.label)

#calculating pathlengths
pathLen=cophenetic.phylo(x = as.phylo(tree2))

#making collecter for max-distance values
#running through each GCF and fishing out internal distances
maxes=data.frame(GCF=unique(df$GCF),maxes= -1)
counter=1
for(i in (unique(df$GCF))) {
  V_indx=which(match(df$GCF,i, )==1)
  Vs=df[V_indx,]

  groupPathLen=pathLen[V_indx,V_indx]
  groupMean=max(groupPathLen[upper.tri(groupPathLen)])

  maxes$maxes[counter]=groupMean
  counter=counter+1
}


maxes=maxes[order(maxes$maxes,decreasing = T),]
df$hiVar=0

pal=c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
"#FFFF99", "#B15928")
for(j in 1:12) {
  df$hiVar[which(match(df$GCF,maxes$GCF[j] )==1)]=j
}
for(k in 13:NROW(maxes)) {
  tree2$tip.label[which(match(df$GCF,maxes$GCF[j] )==1)]=""
}



#plot.phylo(root(tree2,which(pathLen==max(pathLen),arr.ind = T)[1]), tip.color = df$hiVar, type = "fan",show.tip.label = T)

#write.table(x = df,file = "c:/users/milst/Desktop/PseudoforCLC.csv", sep = ";", row.names = F)
write.table(x = df,file = outPath, sep = ";", row.names = F)

