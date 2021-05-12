# RibDif
Count and compare 16S rRNA genes within a genus or species

Takes an input genus and returns a large set of analysed files, the main being phylogenetic trees of full length 16S rRNA genes and common amplicons thereof.

This can be useful for determining if species resolution is possible within certain genera using common metataxonomic approaches


## Installation

#Clone from github

git clone https://github.com/mikaells/RibDif/

#Create conda environment from yml file

conda env create -f RibDif/RibDif.yml

#Activate conda environment 

conda activate RibDif

#Set permissions of files to be executable 

chmod 755 RibDif/*

#silence the parallel citation message if need be (and do with the message what you wish)

parallel --citation   


Note: in_silico_PCR.pl is borrowed from https://github.com/egonozer/in_silico_pcr


## Usage
#Run analysis by as follows 

path-to-RibDif/RibDif.sh -g $genus

#An example run with Ruegeria as target genus:

~/RibDif/RibDif.sh -g Ruegeria

#Can run with species as well, make sure to use quotes

~/rRNA_counter/rRNA_counter.sh -g "Mycoplasma bovis"


Running Ruegeria will take a minute or two, while large genera like Bacillus might take more than an hour as it has 7000+ genomes with 10+ 16S genes.


## Output

The program generates a new directory named after the genus in question. Within that is a summary file and three sub directories.

The main files of interest for amplicon metataxonomic analysis is the amplicon/ folder. The heatmap pdf is the most important one, as it shows which species overlap, which is also summarized in the 'XX-overlap-summary.txt' file.


amplicons/ #files for amplicons
  
    $genus-$region.16S # all 16S genes  
    $genus-$region.aln # alligned from 16S genes (alligned by mafft)  
    $genus-$region.tree # tree from allignment  
    $genus-$region-meta.csv # metadata for each sequence
    $genus-$region-heatmap.pdf # heatmap of 1) clusters found in each genome and 2) genomes that would overlap with given primers
    $region-clusters/ #amplicons clustered into unique clusters by vsearch, in fasta format
    
full/ # files for full length 16S analysis
  
    $genus.16S # all 16S genes  
    $genus.aln # alligned from 16S genes (alligned by mafft)  
    $genus.tree # tree from allignment  
    $genus-meta.csv # metadata for each sequence  

refseq/bacteria/ #contains a folder for each genome, named after refseq-identifier

    GCF*.fna # the original genome file, headers are renamed for downstream compatibility
    GCF*.fna.fai # indexed fna file  
    GCF*.rRNA # rRNA file generated by barrnap  
    GCF*.16S # only the 16S file from .rRNA file  
    GCF*.16sAln # alligned .16S file (by muscle)  
    GCF*.16sTree # tree from .16sAln file (by fasttree)  
    16S_div.pdf # pdf for 16S tree and diversity across genes  
    indiv_16S_dir/ # folder for individual 16S files, for pyani compatibility  
    x.fna #16s files  
    ani/ # files generated from pyani


$genus-summary.csv # contains summary statistics for each genome. Columns 5-8 are in nucleotide mismatches between the genomes' 16S genes, TotalDiv is the cummulative shannon index
 
    #GCF     Genus   Species #16S    Mean    SD      Min     Max   TotalDiv 


  
The tree files are just newick-files and can be used in any treeviewer or R. The meta.csv file is formatted for CLC, but is just a table and can be used for other treeviewers supporting annotation.
