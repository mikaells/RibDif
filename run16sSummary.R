#libraries and plot settings

library(seqinr, quietly = T, warn.conflicts = FALSE)
library(zoo, quietly = T, warn.conflicts = FALSE)
library(ape, quietly = T, warn.conflicts = FALSE)

args=commandArgs(trailingOnly = TRUE);
par(mar = c(2.5, 2.5, 1.8,.5), family="serif", mfrow=c(2,1),mgp = c(1.3, 0.3, 0), font.lab=2)


mismatchPath  =args[1]
alignmentPath =args[2]
pngOut        =args[3] 
fastaPath     =args[4]
treeFile      =args[5]


fullname  = readLines(alignmentPath,n=1)
splitname = strsplit(fullname, split='_')[[1]]

GCF=gsub(">","",paste(splitname[1:2], collapse='_'))
NZ=paste(splitname[3:4], collapse='_')
gen=splitname[5]
spec=paste(splitname[6:length(splitname)-1],collapse=' ')



if(file.exists(mismatchPath)) {
  
  
  mismatch=read.table(mismatchPath, skip=1)[,-1]
  n=NROW(mismatch)
  
  if(n>1) {
    
    upperMismatch=mismatch[upper.tri(mismatch)]
    
    mu=round(mean(upperMismatch),2)
    sd=round(sd(  upperMismatch),2)
    
    mx=max(upperMismatch)
    mi=min(upperMismatch)
    
    
    #address of fasta file
    fasta=read.fasta(alignmentPath)
    
    #turning into a matrix
    DNA_mat=do.call(rbind,fasta)
    
    #finding dimensions of data
    max_len=NCOL(DNA_mat)
    n_seqs=NROW(DNA_mat)
    
    #making empty vector for diversity across string
    divs=rep(-1, max_len) 
    
    #looping across all nucleotide positions
    for( i in 1:max_len) {
      #turn each position into a table with of counts for each nucleotide, while including 0s.
      per_pos=factor(DNA_mat[,i], levels=c("A","G", "C", "T", "-"))
      probs=table(per_pos)/n_seqs
      divs[i]=-sum(probs*ifelse(log(probs)==-Inf, 0, log(probs))) #shannon index
    }
    
    totalDiv=sum(divs)
    
    cat(paste(GCF,gen,spec, n, mu,sd,mi,mx, totalDiv,'\n',sep='\t')) 
    
    #calculation rolling mean of diversity
    roll_means_30=rollmean(divs, k = 30)
    
    tree=read.tree(treeFile)
    
    
    pdf(file = pngOut,width = 7.4,height = 5)
    par(mfrow=c(1,2))
    
    plot(divs, pch=16, cex=0.1, type="p", main="Diversity", xlab="Nucleotide position", ylab="Diversity")
    lines(roll_means_30, col=3)
    legend("topright", legend = "Rolling mean", col=3, lty=1)
    plot(tree)
    dummy=suppressMessages(dev.off())
    
  } else {
    cat(paste(GCF,gen,spec, n, 0,0,0,0, 0,'\n',sep='\t')) 
  } 
  
} else {
  cat(paste(GCF,gen,spec, 0, 0,0,0,0, 0,'\n',sep='\t')) 
}
