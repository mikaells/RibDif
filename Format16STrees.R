rm(list=ls())
suppressMessages(library(ape, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(stringr, quietly = T, warn.conflicts = FALSE))
suppressMessages(library(adephylo, quietly = T, warn.conflicts = FALSE))

args=commandArgs(trailingOnly = TRUE);

treePath = args[1]
outPath  = args[2]
ucPath   = args[3]


#read tree
tree=read.tree(treePath)
ucFile=read.table(ucPath,header=T)

#read cluster file
ucFile=read.table(ucPath)

#removing cluster summaries
ucFile=ucFile[-grep("C",ucFile$V1),]

#get tips for getting genes later
tips=tree$tip.label

#fixing tip names
tips1=str_split_fixed(gsub(":","",gsub("_complete_.*","",gsub("_chromosome_.*","",tips))),pattern = "_", 12)
tips1[,5]=gsub("^(.).*","\\1.",tips1[,5])

#Creating metadata for tips
Name    = tips
GCF     = apply( tips1[ , 1:2 ] , 1 , paste , collapse = "_" )
NZ      = apply( tips1[ , 3:4 ] , 1 , paste , collapse = "_" )
Strain  = trimws(apply( tips1[ , 5:12 ] , 1 , paste , collapse = " " ), which = "right")
Species = apply( tips1[ , 5:6 ] , 1 , paste , collapse = " " )
nGene   = gsub(".*_([0-9]{1,2}$)","\\1",tips)
Cluster = ucFile$V2[match(  Name,ucFile$V9)]


#Making the data frame
df=data.frame(Name, GCF, NZ, Strain, Species, nGene,Cluster)

#write to file
write.table(x = df,file = outPath, sep = ";", row.names = F)

